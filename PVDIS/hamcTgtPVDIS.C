//  hamcTgtPVDIS   -- PVDIS   Target
//  D. Wang  R. Michaels    Jan. 2009

#include "hamcTgtPVDIS.h"
#include "hamcExpt.h"
#include "hamcInout.h"
#include "Rtypes.h"
#include "TRandom.h"
#include <string>
#include <vector>

using namespace std;


#ifndef NODICT
ClassImp(hamcTgtPVDIS)
#endif

hamcTgtPVDIS::hamcTgtPVDIS() : hamcTarget("PVDIS")
{
  did_init = kFALSE;
}

hamcTgtPVDIS::~hamcTgtPVDIS() {
}

Int_t hamcTgtPVDIS::Init(hamcExpt *expt) {

  if (did_init) return OK;

  components.push_back(new hamcTgtSlab(
     "aluminum", 0, 0.005, 0.005, 0.089, 27, 13, 25.3, 2.7));
  components.push_back(new hamcTgtSlab(
     "Liquid Deuterium", 1, 0.25, 0.25, 7.45, 2, 1, 1.876, 0.169));
  components.push_back(new hamcTgtSlab(
     "aluminum", 2, 0.005, 0.005, 0.089, 27, 13, 25.3, 2.7));

  expt->inout->AddToNtuple("zscat",&zscatt);
 
  hamcTarget::Setup();
  hamcTarget::Print();


  did_init = kTRUE;

  return OK;
}

