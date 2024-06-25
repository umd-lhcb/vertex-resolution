// Author: Alex Fernez (as always, borrowing things from Yipeng)
//
// Description: Get weights for +/- variations of fit variables (before any smearing!), where variations
//              are determined by a scaling 1+alpha(log(abs(thetaB_reco-thetaB_true))) weighting up/down
//              events with larger abs(thetaB_reco-thetaB_true)

#include <TFile.h>
#include <TH3D.h>
#include <TMath.h>
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

using namespace std;
using ROOT::RDataFrame;
using ROOT::RDF::RNode;

///////////////////
// Configuration //
///////////////////

static const string DST_TEST_BR = "dst_PX";
static const string D0_TEST_BR  = "d0_PX";

static const string B0_BR_PREFIX  = "b0";
static const string B_BR_PREFIX   = "b";

/////////////
// Helpers //
/////////////

vector<string> setBrPrefix(const string prefix, const vector<string> vars) {
  vector<string> result{};
  for (const auto& v : vars) result.emplace_back(prefix + "_" + v);
  return result;
}

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("GetScaledVariationWeights",
                           "get var weights for scaling up/down large delta_thetaB events");
  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    // I/O
    ("i,input", "specify input ntuple", cxxopts::value<string>())
    ("o,output", "specify output ntuple", cxxopts::value<string>())
    ("t,trees", "specify tree names",
     cxxopts::value<vector<string>>()
     ->default_value("TupleB0/DecayTree,TupleBminus/DecayTree"));
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // get options
  auto ntpNameIn  = parsedArgs["input"].as<string>();
  auto ntpNameOut = parsedArgs["output"].as<string>();

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

    if (treeTest == nullptr) {
      cout << "FATAL: Tree " << t << " doesn't exist!" << endl;
      exit(1);
    }

    if (branchExists(treeTest, DST_TEST_BR)) {
      bMeson = B0_BR_PREFIX;
    } else if (branchExists(treeTest, D0_TEST_BR)) {
      bMeson = B_BR_PREFIX;
    } else {
      cout << "No known branch found for D0 nor D*. Exit now..." << endl;
      exit(1);
    }

    // build a dataframe from input ntuple
    auto           df = static_cast<RNode>(RDataFrame(t, ntpNameIn));
    vector<string> outputBrNames{"runNumber", "eventNumber"};

    // compute variation weights
    // 0.074 is from what Greg found; I mostly just use it because it introduces the right amount of variation for
    // our fit to determine what it wants from the scaling
    // scaling variation weights throw the normalization way off, but gets fixed later
    df = df.Define(
        "thetaB_reco", getBTheta,
        setBrPrefix(bMeson, {"ENDVERTEX_X", "OWNPV_X", "ENDVERTEX_Y", "OWNPV_Y",
                             "ENDVERTEX_Z", "OWNPV_Z"}));
    df = df.Define(
        "thetaB_true", getBTrueTheta,
        setBrPrefix(bMeson, {"TRUEP_X", "TRUEP_Y", "TRUEP_Z"}));
    df = df.Define("wvtx_scale_m", "1 - 0.074*log(TMath::Abs(thetaB_reco-thetaB_true))");
    df = df.Define("wvtx_scale_p", "1 + 0.074*log(TMath::Abs(thetaB_reco-thetaB_true))");
    outputBrNames.emplace_back("wvtx_scale_m");
    outputBrNames.emplace_back("wvtx_scale_p");

    cout << "Writing to " << ntpNameOut << endl;
    df.Snapshot(t, ntpNameOut, outputBrNames, writeOpts);
  }
}
