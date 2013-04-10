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

  Int_t islead=1;
  Int_t isleadonly=0;
  Int_t isthinc12=0;

  sdata = expt->inout->GetStrVect("PREX_target");
  if (sdata.size()>=1) {
    strin = sdata[0];
    if (strin.CmpNoCase("thinc12")==0) {
      islead = 0;
      isthinc12 = 1;
    }
    if (strin.CmpNoCase("leadonly")==0) {
      isleadonly = 1;
      isthinc12 = 0;
    }
  }

  cout << "PREX target choices.  islead = "<<islead<<"   isleadonly = "<<isleadonly<<"  isthinc12 = "<<isthinc12<<endl;

  if (isthinc12) {
       cout << "hamcTgtPREX:: Using the thin C12 target only !"<<endl;
       components.push_back(new hamcTgtSlab(
           "diamond", 0, 0.0005, 0.00005, 0.0188, 12, 6, 11.25, 3.52));
  }

  
  sdata = expt->inout->GetStrVect("PREX_model");
  if (sdata.size()>=1) {
    strin = sdata[0];
    if (strin.CmpNoCase("horca40")==0) {
        cout << "Using a Calcium 40 Target "<<endl;
        // Assume RL = 20 g/cm^2 = 12.9 cm.  5% RL = 0.65 cm.  
        // for 5% RL, Tgt eff = 0.37 * tgt (instead of 0.27*) = 0.24 cm
        components.push_back(new hamcTgtSlab(
            // name   id   len    efflen  rlen     a    z    m    dens
            "calcium40", 0, 0.0065, 0.0024, 0.129, 40,  20, 37.5, 1.55));
        islead=0;
    }
    if (strin.CmpNoCase("horca48")==0) {
        cout << "Using a Calcium 48 Target "<<endl;
        // Assume RL = 20 g/cm^2 = 12.9 cm.  5% RL = 0.65 cm.  
        // for 5% RL, Tgt eff = 0.37 * tgt (instead of 0.27*) = 0.24 cm
        components.push_back(new hamcTgtSlab(
            // name   id   len    efflen  rlen   a  z    m   dens
            "calcium48", 0, 0.0065, 0.0024, 0.129, 48,  20, 45, 1.55));
       islead=0;
    }
    if (strin.CmpNoCase("horsn")==0) {
       cout << "Using a Tin Target "<<endl;
       // need to check RL and density
       if (!isthinc12) components.push_back(new hamcTgtSlab(
          "tin", 0, 0.0016, 0.00059, 1.6, 120, 50, 112, 5.8));
     islead=0;
    }
  }

  // (updated Sept 17, was 0.37, now 0.27)
  // Tgt_eff = 0.27 * Tgt_Len, accounts for radiative tail
  // not already in hamc event generation.

  if (islead && !isthinc12) {
    if (!isleadonly) components.push_back(new hamcTgtSlab(
       // name  id   len        elen     rlen   a   z   m      d
      "diamond", 0, 0.00015, 0.000041, 0.0188, 12, 6, 11.25, 3.52));

    components.push_back(new hamcTgtSlab(
      "lead",    1, 0.00050, 0.000135, 0.0056, 208, 82, 195, 11.35));

    if (!isleadonly) components.push_back(new hamcTgtSlab(
      "diamond", 0, 0.00015, 0.000041, 0.0188, 12, 6, 11.25, 3.52));

  }

  expt->inout->AddToNtuple("zscat",&zscatt);

  hamcTarget::Setup();
  hamcTarget::Print();

  did_init = kTRUE;

  return OK;
}

