//  hamcTgtPREX   -- PREX  Pb/C   Target
//  R. Michaels  June 2008

#include "hamcTgtPREX.h"
#include "hamcExpt.h"
#include "hamcInout.h"
#include "Rtypes.h"
#include "TRandom.h"
#include <string>
#include <vector>

using namespace std;


#ifndef NODICT
ClassImp(hamcTgtPREX)
#endif

hamcTgtPREX::hamcTgtPREX() : hamcTarget("PREX")
{
  did_init = kFALSE;
}

hamcTgtPREX::~hamcTgtPREX() {
}

Int_t hamcTgtPREX::Init(hamcExpt *expt) {

  if (did_init) return OK;

  components.push_back(new hamcTgtSlab(
     "diamond", 0, 0.015, 0.005, 18.8, 12, 6, 11.25, 3.52));
  components.push_back(new hamcTgtSlab(
     "lead",    1, 0.050, 0.017, 0.56, 208, 82, 195, 11.35));
  components.push_back(new hamcTgtSlab(
     "diamond", 0, 0.015, 0.005, 18.8, 12, 6, 11.25, 3.52));

  expt->inout->AddToNtuple("zscat",&zscatt);

  hamcTarget::Setup();
  hamcTarget::Print();

  did_init = kTRUE;

  return OK;
}

