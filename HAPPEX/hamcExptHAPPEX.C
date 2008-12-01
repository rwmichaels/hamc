//  hamcExptHAPPEX   -- HAPPEX Experiment.
//  R. Michaels  Dec 2008

#include "hamcExptHAPPEX.h"
#include "hamcExpt.h"
#include "hamcSpecHRS.h"
#include "hamcTgtHAPPEX.h"
#include "hamcPhyHAPPEX.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptHAPPEX)
#endif

using namespace std;

hamcExptHAPPEX::hamcExptHAPPEX() : hamcSingles("HAPPEX")
{
  physics = new hamcPhyHAPPEX();
  target  = new hamcTgtHAPPEX();
}

hamcExptHAPPEX::~hamcExptHAPPEX() {
}

Int_t hamcExptHAPPEX::Init(string sfile) {

// Defaults in case the setup file is empty
  hamcSingles::SetP0(3.176);
  hamcSingles::SetTheta(6.106);  // for R-HRS

  hamcSingles::Init(sfile);

  // Defaults (can be over-ridden by input file)
  static_cast<hamcSpecHRS*>(spectrom[0])->UseCollimator();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseColdSeptum();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseMatrixTrans();

  spectrom[0]->Init(this);
  spectrom[0]->Print();

  return OK;
}

void hamcExptHAPPEX::Analysis() {

// Fill some histograms for diagnostic purposes.

  Int_t lprint = 1;  // to turn on(1) or off(0) local print

  if (lprint) cout << "Into hamcExptHAPPEX:: Analysis "<<endl;

}

