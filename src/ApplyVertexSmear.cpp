// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Nov 21, 2022 at 09:57 PM -0500
//
// Description: Apply vertex smearing to ntuples

// #include <Math/Vector3D.h>
// #include <Math/Vector4D.h>
// #include <Math/VectorUtil.h>
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

static const string B0_BR_PREFIX  = "b0";
static const string B_BR_PREFIX   = "b";
static const string DST_BR_PREFIX = "dst";
static const string D0_BR_PREFIX  = "d0";

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

/////////////////////////////////
// Rest frame variable helpers //
/////////////////////////////////

vector<float> loadDeltaTheta(string auxFile) {
  vector<float> result{};
  auto          df = RDataFrame(THETA_TREE_NAME, auxFile);
  df.Foreach(
      [&](float x) {
        // if (x > 0.25 || x < -0.25) return;
        // result.emplace_back(TMath::Abs(x));
        result.emplace_back(x);
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

auto computeDeltaThetaHelper(double& lin, double& quad) {
  auto rng = make_shared<TRandomMixMax256>(RAND_SEED * 2);

  return [&lin, &quad, rng](float rawAngle) {
    rawAngle = TMath::Abs(rawAngle);
    int sign = 1;
    if (rng->Uniform(0, 1) > 0.5) sign = -1;
    return static_cast<float>(sign *
                              (lin * rawAngle + quad * rawAngle * rawAngle));
  };
}

vector<string> setBrPrefix(const string prefix, const vector<string> vars) {
  vector<string> result{};
  for (const auto& v : vars) result.emplace_back(prefix + "_" + v);
  return result;
}

RNode computeFitVars(RNode df, double mB, string bMeson, string dMeson,
                     string suffix, vector<string>& outputBrs) {
  auto estB4MomPartial = [=](PxPyPzEVector v4BReco, XYZVector v3BFlight) {
    return estB4Mom(v4BReco, v3BFlight, mB);
  };
  suffix = "_" + suffix;

  vector<string> brNames = {"mm2", "q2", "el"};
  for (auto& n : brNames) outputBrs.emplace_back(n + suffix);

  auto fourVecHelper = [](double px, double py, double pz, double pe) {
    return PxPyPzEVector(px, py, pz, pe);
  };
  vector<string> kinBrNames = {"PX", "PY", "PZ", "PE"};

  return df.Define("v4_b_reco", fourVecHelper, setBrPrefix(bMeson, kinBrNames))
      .Define("v4_d", fourVecHelper, setBrPrefix(dMeson, kinBrNames))
      .Define("v4_mu", fourVecHelper, setBrPrefix("mu", kinBrNames))
      .Define("v4_b_est", estB4MomPartial, {"v4_b_reco", "v3_b_dir"})
      .Define("mm2" + suffix, m2Miss, {"v4_b_est", "v4_b_reco"})
      .Define("q2" + suffix, q2, {"v4_b_est", "v4_d"})
      .Define("el" + suffix, el, {"v4_b_est", "v4_mu"});
}

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

    // find B meson branch name
    auto   ntpTest  = TFile(ntpNameIn.c_str());
    auto   treeTest = static_cast<TTree*>(ntpTest.Get(t.c_str()));
    string bMeson;
    string dMeson;
    double mB;

    if (treeTest == nullptr) {
      cout << "FATAL: Tree " << t << " doesn't exist!" << endl;
      exit(1);
    }

    if (branchExists(treeTest, DST_TEST_BR)) {
      bMeson = B0_BR_PREFIX;
      dMeson = DST_BR_PREFIX;
      mB     = B0_M;
    } else if (branchExists(treeTest, D0_TEST_BR)) {
      bMeson = B_BR_PREFIX;
      dMeson = D0_BR_PREFIX;
      mB     = B_M;
    } else {
      cout << "No known branch found for D0 nor D*. Exit now..." << endl;
      exit(1);
    }

    // reinitialize random seed for each tree
    auto funcSmr = getRandSmrHelper(vDeltaTheta);

    // build a dataframe from input ntuple
    auto           df = static_cast<RNode>(RDataFrame(t, ntpNameIn));
    vector<string> outputBrNames{"runNumber", "eventNumber"};

    // define raw branches
    df = defineBranch(df, "");
    for (auto& [br, expr] : FIT_VARS) outputBrNames.emplace_back(br);

    // get a random angle
    df = df.Define("raw_delta_theta", funcSmr, {});
    outputBrNames.emplace_back("raw_delta_theta");

    auto funcAngle = computeDeltaThetaHelper(fitLin, fitQuad);
    df = df.Define(bMeson + "_delta_theta", funcAngle, {"raw_delta_theta"});
    outputBrNames.emplace_back(bMeson + "_delta_theta");

    // smear B meson flight vector angle theta
    df = df.Define(
        "v3_b_dir", buildBFlightDir,
        setBrPrefix(bMeson, {"ENDVERTEX_X", "OWNPV_X", "ENDVERTEX_Y", "OWNPV_Y",
                             "ENDVERTEX_Z", "OWNPV_Z", "delta_theta"}));

    // recompute fit vars
    df = computeFitVars(df, mB, bMeson, dMeson, "vtx_smr", outputBrNames);

    // compute optional variation weights
    // the signs ARE correct
    df = df.Define("wvtx_debug",
                   "0.01*log(TMath::Abs(" + bMeson + "_delta_theta))");
    df = df.Define("wvtx_m", "1 + wvtx_debug");
    df = df.Define("wvtx_p", "1 - wvtx_debug");
    outputBrNames.emplace_back("wvtx_debug");
    outputBrNames.emplace_back("wvtx_p");
    outputBrNames.emplace_back("wvtx_m");

    cout << "Writing to " << ntpNameOut << endl;
    df.Snapshot(t, ntpNameOut, outputBrNames, writeOpts);
  }
}
