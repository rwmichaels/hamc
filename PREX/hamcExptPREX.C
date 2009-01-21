//  hamcExptPREX   -- Pb Radius Experiment.
//  R. Michaels  June 2008

#include "hamcExptPREX.h"
#include "hamcExpt.h"
#include "hamcSingles.h"
#include "hamcSpecHRS.h"
#include "hamcTgtPREX.h"
#include "hamcPhyPREX.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptPREX)
#endif

using namespace std;

hamcExptPREX::hamcExptPREX() : hamcSingles("PREX")
{
  physics = new hamcPhyPREX();
  target  = new hamcTgtPREX();
}

hamcExptPREX::~hamcExptPREX() {
}

Int_t hamcExptPREX::Init(string sfile) {

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

