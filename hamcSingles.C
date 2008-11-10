//  hamcSingles   -- base class for single-arm experiment
//  R. Michaels  Sept 2007

#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcSpecHRS.h"
#include "hamcInout.h"
#include "THaString.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcSingles)
#endif


hamcSingles::hamcSingles(string sname) : hamcExpt(sname)
{
  event = new hamcEvent();  // should perhaps be a SinglesEvent ?
}

hamcSingles::~hamcSingles() {
   if (event) delete event;
}


Int_t hamcSingles::Init(string sfile) {

  hamcExpt::InitInput(sfile);
  THaString strin;
  Int_t which = RIGHTHRS; // default

  vector<string> sdata; 
  sdata = inout->GetStrVect("HRS_arm");
  if (sdata.size()>=1) {
    strin = sdata[0];
    if (strin.CmpNoCase("left")==0) {
      cout << "HRS arm = left "<<endl;
      which = LEFTHRS;
    }
    if (strin.CmpNoCase("right")==0) {
      cout << "HRS arm = right "<<endl;
      which = RIGHTHRS;
    }
  }
  sdata = inout->GetStrVect("HRS_P0");
  if (sdata.size()>=1) {
     sscanf(sdata[0].c_str(),"%f",&P0);
  }
  sdata = inout->GetStrVect("HRS_angle");
  if (sdata.size()>=1) {
     sscanf(sdata[0].c_str(),"%f",&angle);
  }

  cout << "hamcSingles:  P0 = "<<P0<<"    angle = "<<angle<<endl;

  SetSpectrom(which, P0, angle);

  return hamcExpt::Init(sfile);

}

Int_t hamcSingles::SetSpectrom(Int_t which, Float_t pmom, Float_t theta) {

  spectrom.push_back(new hamcSpecHRS(which, pmom, theta));

  return OK;
}




   
