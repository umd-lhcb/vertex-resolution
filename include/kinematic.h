// Author: Yipeng Sun, Svende Braun
// License: BSD 2-clause
// Last Change: Mon Nov 21, 2022 at 12:21 PM -0500

#pragma once

#include <vector>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TROOT.h>
#include <TVector3.h>

using std::vector;

using ROOT::Math::PxPyPzEVector;
using ROOT::Math::PxPyPzMVector;
using ROOT::Math::XYZVector;

#define K_M 493.677
#define PI_M 139.570
#define B_M 5279.34
#define B0_M 5279.65

//////////////////////
// Rebuild momentum //
//////////////////////

XYZVector buildBFlightDir(double endVtxX, double ownPvX, double endVtxY,
                          double ownPvY, double endVtxZ, double ownPvZ,
                          float smrAngle) {
  // FIXME: Legacy way of doing things
  auto flight = TVector3(endVtxX - ownPvX, endVtxY - ownPvY, endVtxZ - ownPvZ);
  flight.SetTheta(flight.Theta() + smrAngle);
  return XYZVector(flight);
}

//////////////////////////////
// Rest frame approximation //
//////////////////////////////

PxPyPzEVector estB4Mom(PxPyPzEVector v4BReco, XYZVector v3BFlight,
                       double mBRef = B_M) {
  auto mB  = v4BReco.M();
  auto pzB = v4BReco.Pz();

  auto cosX = v3BFlight.Unit().X();
  auto cosY = v3BFlight.Unit().Y();
  auto cosZ = v3BFlight.Unit().Z();

  Double_t pBMag = (mBRef / mB) * pzB / cosZ;
  return PxPyPzEVector(pBMag * cosX, pBMag * cosY, pBMag * cosZ,
                       TMath::Sqrt(pBMag * pBMag + mBRef * mBRef));
}

// all in GeV(^2)!
// also removed all template parameters because RDataFrame doesn't like them.
Double_t m2Miss(PxPyPzEVector& v4BEst, PxPyPzEVector& v4BReco) {
  return (v4BEst - v4BReco).M2() / 1000 / 1000;
}

Double_t el(PxPyPzEVector& v4BEst, PxPyPzEVector& v4Mu) {
  auto boost    = v4BEst.BoostToCM();
  auto v4MuRest = ROOT::Math::VectorUtil::boost(v4Mu, boost);
  return v4MuRest.E() / 1000;
}

Double_t q2(PxPyPzEVector& v4BEst, PxPyPzEVector& v4D) {
  return (v4BEst - v4D).M2() / 1000 / 1000;
}

// in MeV!
Double_t calcBM(PxPyPzEVector& v4BReco) { return v4BReco.M(); }
