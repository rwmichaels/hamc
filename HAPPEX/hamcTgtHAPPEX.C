//  hamcTgtHAPPEX   -- HAPPEX   Target
//  R. Michaels  June 2008

#include "hamcTgtHAPPEX.h"
#include "hamcExpt.h"
#include "hamcInout.h"
#include "Rtypes.h"
#include "TRandom.h"
#include <string>
#include <vector>

using namespace std;


#ifndef NODICT
ClassImp(hamcTgtHAPPEX)
#endif

hamcTgtHAPPEX::hamcTgtHAPPEX() : hamcTarget("HAPPEX")
{
  did_init = kFALSE;
}

hamcTgtHAPPEX::~hamcTgtHAPPEX() {
}

Int_t hamcTgtHAPPEX::Init(hamcExpt *expt) {

  if (did_init) return OK;

  components.push_back(new hamcTgtSlab(
     "aluminum", 0, 0.001, 0.0014, 0.089, 14, 13, 25.3, 2.7));
  components.push_back(new hamcTgtSlab(
     "hydrogen", 1, 0.00050, 0.2, 8.66, 1, 0, 0.938, 0.0708));
  components.push_back(new hamcTgtSlab(
     "aluminum", 2, 0.001, 0.00178, 0.089, 14, 13, 25.3, 2.7));

  // FIXME : may move this to base class
  expt->inout->AddToNtuple("zscat",&zscatt);

  hamcTarget::Setup();
  hamcTarget::Print();

  did_init = kTRUE;

  return OK;
}

