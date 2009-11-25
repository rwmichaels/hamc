//  hamcTgtPREX   -- PREX  Pb/C   Target
//  R. Michaels  June 2008

#include "hamcTgtPREX.h"
#include "hamcExpt.h"
#include "hamcInout.h"
#include "THaString.h"
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

  THaString strin;
  vector<string> sdata; 
  sdata = expt->inout->GetStrVect("PREX_model");
  Int_t islead=1;
  if (sdata.size()>=1) {
    strin = sdata[0];
    if (strin.CmpNoCase("horca")==0) {
      cout << "Using a Calcium Target "<<endl;
      // need to check RL and density
     components.push_back(new hamcTgtSlab(
	  // 10% of RL = 0.66; Teff = 0.37*T (RL loss)
       "calcium", 0, 0.0066, 0.0024, 6.6, 48, 20, 45, 1.6));
     islead=0;
    }
    if (strin.CmpNoCase("horsn")==0) {
      cout << "Using a Tin Target "<<endl;
      // need to check RL and density
     components.push_back(new hamcTgtSlab(
       "tin", 0, 0.0016, 0.00059, 1.6, 120, 50, 112, 5.8));
     islead=0;
    }
  }

  // Tgt_eff = 0.37 * Tgt_Len, accounts for radiative tail
  // not already in hamc event generation.
  
  if (islead) {
    components.push_back(new hamcTgtSlab(
      "diamond", 0, 0.00015, 0.000056, 0.0188, 12, 6, 11.25, 3.52));
    components.push_back(new hamcTgtSlab(
      "lead",    1, 0.00050, 0.000186, 0.0056, 208, 82, 195, 11.35));
    components.push_back(new hamcTgtSlab(
      "diamond", 0, 0.00015, 0.000056, 0.0188, 12, 6, 11.25, 3.52));
  }

  expt->inout->AddToNtuple("zscat",&zscatt);

  hamcTarget::Setup();
  hamcTarget::Print();

  did_init = kTRUE;

  return OK;
}

