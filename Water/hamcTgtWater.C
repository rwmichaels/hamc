//  hamcTgtWater   -- Water Cell Target
//  R. Michaels  June 2008

#include "hamcTgtWater.h"
#include "hamcExpt.h"
#include "hamcInout.h"
#include "Rtypes.h"
#include "TRandom.h"
#include <string>
#include <vector>

using namespace std;


#ifndef NODICT
ClassImp(hamcTgtWater)
#endif

hamcTgtWater::hamcTgtWater() : hamcTarget("Water")
{
  did_init = kFALSE;
}

hamcTgtWater::~hamcTgtWater() {
}

Int_t hamcTgtWater::Init(hamcExpt *expt) {

  if (did_init) return OK;

  components.push_back(new hamcTgtSlab(
       // name   id   len    elen   rlen   a   z   m      d
      "oxygen",   0,  0.005,   0.0045,  36.1, 16,  8, 15.0, 0.88));
   components.push_back(new hamcTgtSlab(
      "hydrogen", 1,  0.005,   0.0045,  36.1,  1,  1, 0.938, 0.12)); 

  expt->inout->AddToNtuple("zscat",&zscatt);
 
  hamcTarget::Setup();
  hamcTarget::Print();


  did_init = kTRUE;

  return OK;
}

