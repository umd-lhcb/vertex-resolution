// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Nov 21, 2022 at 06:11 AM -0500
//
// Description: Apply vertex smearing to ntuples

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TFile.h>
#include <TH3D.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandomGen.h>
#include <TString.h>
#include <TTree.h>

#include <ROOT/RDataFrame.hxx>
#include <boost/range/join.hpp>
#include <cxxopts.hpp>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "kinematic.h"
#include "utils.h"

#define RAND_SEED 42

using namespace std;
using ROOT::RDataFrame;
using ROOT::Math::PxPyPzEVector;
using ROOT::Math::XYZVector;
using ROOT::RDF::RNode;

///////////////////
// Configuration //
///////////////////

typedef vector<pair<string, string>> vPStrStr;

static const vPStrStr FIT_VARS{
    {"q2_input", "FitVar_q2 / 1000 / 1000"},
    {"mm2_input", "FitVar_Mmiss2 / 1000 / 1000"},
    {"el_input", "FitVar_El / 1000"},
};

static const string DST_TEST_BR = "dst_PX";
static const string D0_TEST_BR  = "d0_PX";

static const string B0_BR_PREFIX = "b0";
static const string B_BR_PREFIX  = "b";

static const string THETA_TREE_NAME = "Smear";
static const string THETA_BR_NAME   = "Delta";

////////////////////////
// RDataFrame helpers //
////////////////////////

// Idea stolen from:
//   https://root-forum.cern.ch/t/running-rdataframes-define-in-for-loop/32484/2
RNode defineBranch(RNode df, string particle = B0_BR_PREFIX,
                   const vPStrStr& rules = FIT_VARS, int idx = 0) {
  if (rules.size() == idx) return df;

  auto inputBrName = rules[idx].second;
  if (particle != ""s) inputBrName = particle + "_" + inputBrName;
  cout << "Define " << rules[idx].first << " as " << inputBrName << endl;

  return defineBranch(df.Define(rules[idx].first, inputBrName), particle, rules,
                      idx + 1);
}

////////////////////////////////////
// Helpers for weight computation //
////////////////////////////////////

// template <typename T, typename C = decay_t<decltype(*begin(declval<T>()))>>
// tuple<RNode, vector<string>, vector<TH3D*>> applyWtFromHistos(
//     RNode df, TFile* ntpHisto, string histoPrefix, string weightBrPrefix,
//     T iterable) {
//   auto outputBrs = vector<string>{};
//   auto histos    = vector<TH3D*>{};

//   for (const auto& h : iterable) {
//     auto histoName = string(histoPrefix + "__" + h);
//     auto histoWt   = static_cast<TH3D*>(ntpHisto->Get(histoName.data()));
//     histos.emplace_back(histoWt);
//     cout << "  Loading histo " << histoName << endl;

//     double prescale = 1.0;

//     auto brName = weightBrPrefix + "_" + h;
//     cout << "  Generating " << brName << "..." << endl;
//     df = df.Define(brName,
//                    [histoWt, prescale](double& x, double& y, double& z) {
//                      auto binIdx = histoWt->FindFixBin(x, y, z);
//                      return histoWt->GetBinContent(binIdx) * prescale;
//                    },
//                    {"P", "ETA", "nTracks"});
//     outputBrs.emplace_back(brName);
//   }

//   return {df, outputBrs, histos};
// }

// pair<vPStrStr, vector<string>> genWtDirective(YAML::Node    node,
//                                               const string& wtPrefix,
//                                               string brPrefix = "is_misid_")
//                                               {
//   vPStrStr       directives{};
//   vector<string> outputBrs{};
//   const auto     wtTargetParticle = "MuTag";
//   const auto     wtSmrParticles   = {"k", "pi"};

//   vector<string> particles{};
//   // first find particles
//   for (auto it = node.begin(); it != node.end(); it++)
//     particles.emplace_back(it->first.as<string>());

//   // generate the automatic weight for each event based on the species of the
//   // event
//   auto expr  = ""s;
//   auto first = true;

//   for (const auto& p : particles) {
//     auto wtBrName = wtPrefix + "_" + p + "TagTo" + wtTargetParticle;
//     if (!first) expr += " + ";
//     first = false;
//     expr += brPrefix + p + "*" + wtBrName;
//   }
//   outputBrs.emplace_back(wtPrefix);
//   directives.emplace_back(pair{wtPrefix, expr});
//   cout << "  " << wtPrefix << " = " << expr << endl;

//   // generate the DiF smearing weight for each event
//   vector<string> brSmrNames{};
//   for (const auto& tgt : wtSmrParticles) {
//     expr  = ""s;
//     first = true;
//     for (const auto& p : particles) {
//       auto wtBrName = wtPrefix + "_" + p + "TagTo" + capitalize(tgt) +
//       "True"; if (!first) expr += " + "; first = false; expr += brPrefix + p
//       + "*" + wtBrName;
//     }

//     auto outputBr = wtPrefix + "_smr_" + tgt;
//     brSmrNames.push_back(outputBr);
//     outputBrs.emplace_back(outputBr);
//     directives.emplace_back(pair{outputBr, expr});
//     cout << "  " << outputBr << " = " << expr << endl;
//   }

//   // generate the DiF no smearing weight
//   auto brNoSmr = wtPrefix + "_no_smr";
//   outputBrs.emplace_back(wtPrefix + "_no_smr");

//   expr = "1.0"s;
//   for (const auto& smr : brSmrNames) {
//     expr += " - " + smr;
//   }
//   cout << "  " << brNoSmr << " = " << expr << endl;
//   directives.emplace_back(pair{brNoSmr, expr});

//   return {directives, outputBrs};
// }

/////////////////////////////////
// Rest frame variable helpers //
/////////////////////////////////

vector<float> loadDeltaTheta(string auxFile) {
  vector<float> result{};
  auto          df = RDataFrame(THETA_TREE_NAME, auxFile);
  df.Foreach(
      [&](float x) {
        if (x > 0.25 || x < -0.25) return;
        result.emplace_back(TMath::Abs(x));
      },
      {THETA_BR_NAME});
  return result;
}

auto getRandSmrHelper(vector<float>& smr) {
  auto size = make_shared<unsigned long>(smr.size());
  auto rng  = make_shared<TRandomMixMax256>(RAND_SEED);

  return [&smr, size, rng] {
    unsigned long rand = rng->Uniform(0, *(size.get()));
    return smr[rand];
  };
}

// template <typename F>
// RNode computeDiFVars(RNode df, F& randGetter, double mB, string suffix,
//                      vector<string>& outputBrs) {
//   // we probably did some unnecessary copies here, but deducing those nested
//   // lambdas can be quite hard so I'm just being lazy here.
//   auto rebuildMu4MomPartial = [=, &randGetter](PxPyPzEVector v4Mu) {
//     vector<double> smr = randGetter();
//     return rebuildMu4Mom(v4Mu, smr);
//   };
//   auto estB4MomPartial = [=](PxPyPzEVector v4BReco, XYZVector v3BFlight) {
//     return estB4Mom(v4BReco, v3BFlight, mB);
//   };

//   vector<string> brNames = {"mm2", "q2", "el", "b_m"};
//   for (auto& n : brNames) outputBrs.emplace_back(n + suffix);

//   return df.Define("v4_mu" + suffix, rebuildMu4MomPartial, {"v4_mu"})
//       .Define("v4_b_reco" + suffix, "v4_mu" + suffix + " + v4_d")
//       .Define("v4_b_est" + suffix, estB4MomPartial,
//               {"v4_b_reco" + suffix, "v3_b_dir"})
//       .Define("mm2" + suffix, m2Miss,
//               {"v4_b_est" + suffix, "v4_b_reco" + suffix})
//       .Define("q2" + suffix, q2, {"v4_b_est" + suffix, "v4_d"})
//       .Define("el" + suffix, el, {"v4_b_est" + suffix, "v4_mu" + suffix})
//       .Define("b_m" + suffix, calcBM, {"v4_b_reco" + suffix});
// }

// template <typename F1, typename F2>
// pair<RNode, vector<string>> defRestFrameVars(RNode df, TTree* tree,
//                                              F1& randKGetter,
//                                              F2& randPiGetter) {
//   vector<string> outputBrs{};
//   string         dMeson = ""s;
//   string         bMeson = ""s;
//   double         mBRef;

//   if (branchExists(tree, DST_TEST_BR)) {
//     bMeson = B0_BR_PREFIX;
//   } else if (branchExists(tree, D0_TEST_BR)) {
//     bMeson = B_BR_PREFIX;
//   } else {
//     cout << "No known branch found for D0 nor D*. Exit now..." << endl;
//     exit(1);
//   }

//   // Basic vectors that we use and not going to change
//   df = df.Define(
//              "v4_d",
//              [](double px, double py, double pz, double e) {
//                return PxPyPzEVector(px, py, pz, e);
//              },
//              setBrPrefix(dMeson, {"PX", "PY", "PZ", "PE"}))
//            .Define(
//                "v4_mu",
//                [](double px, double py, double pz, double e) {
//                  return PxPyPzEVector(px, py, pz, e);
//                },
//                setBrPrefix(MU_BR_PREFIX, {"PX", "PY", "PZ", "PE"}))
//            .Define("v3_b_dir", buildBFlightDir,
//                    setBrPrefix(bMeson, {"ENDVERTEX_X", "OWNPV_X",
//                    "ENDVERTEX_Y",
//                                         "OWNPV_Y", "ENDVERTEX_Z",
//                                         "OWNPV_Z"}));

//   // Replace mass hypo and compute fit vars
//   df = computeDiFVars(df, randPiGetter, mBRef, "_smr_pi", outputBrs);
//   df = computeDiFVars(df, randKGetter, mBRef, "_smr_k", outputBrs);

//   return {df, outputBrs};
// }

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("ApplyMisIDWeight",
                           "unfolding weihgts applyer (A).");

  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    // I/O
    ("i,input", "specify input ntuple", cxxopts::value<string>())
    ("x,aux", "specify auxiliary ntuple", cxxopts::value<string>())
    ("o,output", "specify output ntuple", cxxopts::value<string>())
    ("t,trees", "specify tree names",
     cxxopts::value<vector<string>>()
     ->default_value("TupleB0/DecayTree,TupleBminus/DecayTree"))
    // fit params
    ("fitLin", "specify linear coeff", cxxopts::value<double>()
     ->default_value("0.105"))
    ("fitQuad", "specify quadratic coeff", cxxopts::value<double>()
     ->default_value("6.29"))
  ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // get options
  auto ntpNameIn  = parsedArgs["input"].as<string>();
  auto ntpNameOut = parsedArgs["output"].as<string>();
  auto ntpNameAux = parsedArgs["aux"].as<string>();

  auto fitLin  = parsedArgs["fitLin"].as<double>();
  auto fitQuad = parsedArgs["fitQuad"].as<double>();

  // load true and reco'ed flight theta angles
  auto vDeltaTheta = loadDeltaTheta(ntpNameAux);

  // snapshot option
  auto writeOpts  = ROOT::RDF::RSnapshotOptions{};
  writeOpts.fMode = "UPDATE";

  // loop over input trees
  auto inputTrees = parsedArgs["trees"].as<vector<string>>();
  for (auto& t : inputTrees) {
    cout << "--------" << endl;
    cout << "Working on tree: " << t << endl;

    // reinitialize random seed for each tree
    auto funcSmr = getRandSmrHelper(vDeltaTheta);

    // build a dataframe from input ntuple
    auto           df = static_cast<RNode>(RDataFrame(t, ntpNameIn));
    vector<string> outputBrNames{"runNumber", "eventNumber"};

    // define raw branches
    df = defineBranch(df, "");
    for (auto& [br, expr] : FIT_VARS) outputBrNames.emplace_back(br);

    cout << "Writing to " << ntpNameOut << endl;
    df.Snapshot(t, ntpNameOut, outputBrNames, writeOpts);
  }

  // for (auto it = outputDirective.begin(); it != outputDirective.end(); it++)
  // {
  //   auto treeName    = it->first.as<string>();
  //   auto ntpsToClean = vector<TFile*>{ntpInTest};

  //   if (applyAlias) {
  //     df = defineBranch(df, particle);
  //     // compute ETA
  //     df = df.Define("ETA",
  //                    [](double& p, double& pz) {
  //                      return 0.5 * TMath::Log((p + pz) / (p - pz));
  //                    },
  //                    {"P", "PZ"});
  //   }

  //   // read smearing factors from aux ntuple
  //   auto smrFacK  = vector<vector<double>>{};
  //   auto smrFacPi = vector<vector<double>>{};
  //   getSmrFac(smrFacK, ntpAux, kSmrBrName);
  //   getSmrFac(smrFacPi, ntpAux, piSmrBrName);

  //   // we can call these functions directly to get random smearing factors
  //   auto randSmrFacK  = getRandSmrHelper(smrFacK);
  //   auto randSmrFacPi = getRandSmrHelper(smrFacPi);

  //   // Recompute fit vars
  //   auto [dfOut, outputBrsFitVars] =
  //       defRestFrameVars(df, treeTest, randSmrFacK, randSmrFacPi);
  //   for (auto br : outputBrsFitVars) outputBrNames.emplace_back(br);
  //   df = dfOut;

  //   // add species tags
  //   auto [directivesTags, outputBrsTags] =
  //       genTaggedCutDirective(ymlConfig["tags"]);
  //   df = defineBranch(df, ""s, directivesTags);
  //   for (const auto& br : outputBrsTags) outputBrNames.emplace_back(br);

  //   // add all kinds of weights
  //   for (auto entry : it->second) {
  //     auto histoPrefix    = entry["prefix"].as<string>();
  //     auto histoFile      = entry["file"].as<string>();
  //     auto weightBrPrefix = entry["name"].as<string>();
  //     histoFile           = filePrefix + "/" + histoFile;

  //     cout << "Handling tree " << treeName << " from histos of prefix "
  //          << histoPrefix << " from file " << histoFile << endl;

  //     // add weights required by misID weights
  //     auto ntpHisto = new TFile(histoFile.data(), "READ");
  //     ntpsToClean.emplace_back(ntpHisto);

  //     auto histoWtNames    = buildHistoWtNames(particle, ymlConfig["tags"]);
  //     auto histoSmrWtNames = buildHistoSmrWtnames(ymlConfig["tags"]);
  //     cout << "Generate transfer factors/DiF smearing wieghts for all
  //     species"
  //          << endl;
  //     auto [dfHistos, outputBrsHistos, histos] =
  //         applyWtFromHistos(df, ntpHisto, histoPrefix, weightBrPrefix,
  //                           boost::join(histoWtNames, histoSmrWtNames));
  //     for (auto br : outputBrsHistos) outputBrNames.emplace_back(br);

  //     // add the actual misID weights
  //     auto [directives, outputBrsWts] =
  //         genWtDirective(ymlConfig["tags"], weightBrPrefix);
  //     df = defineBranch(dfHistos, ""s, directives);
  //     for (auto br : outputBrsWts) outputBrNames.emplace_back(br);
  //   }
}
