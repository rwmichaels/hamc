//  hamcExptPVDIS   -- Pb Radius Experiment.
//  D. Wang   R. Michaels  Jan. 2009

#include "hamcExptPVDIS.h"
#include "hamcExpt.h"
#include "hamcSingles.h"
#include "hamcSpecHRS.h"
#include "hamcTgtPVDIS.h"
#include "hamcPhyPVDIS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptPVDIS)
#endif

using namespace std;

hamcExptPVDIS::hamcExptPVDIS() : hamcSingles("PVDIS")
{
  physics = new hamcPhyPVDIS();
  target  = new hamcTgtPVDIS();
}

hamcExptPVDIS::~hamcExptPVDIS() {
}

Int_t hamcExptPVDIS::Init(string sfile) {

// Defaults in case the setup file is empty
  hamcSingles::SetP0(1.04);
  hamcSingles::SetTheta(5.0);

  hamcSingles::Init(sfile);

  // Defaults (can be over-ridden by input file)
  static_cast<hamcSpecHRS*>(spectrom[0])->UseCollimator();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseWarmSeptum();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseMatrixTrans();

  spectrom[0]->Init(this);
  spectrom[0]->Print();

  return OK;
}

