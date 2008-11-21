//  hamcRad   -- Internal and external Brehmstrahlung
//  This will be a member of hamcPhysics.

//  Will get rid of the "NOTSTANDALONE" ifdefs soon.
//  This illustrates how to develop a new class, first 
//  isolation, and then coupled to the other classes.

//  R. Michaels  Nov 2008


#ifdef NOTSTANDALONE
#include "hamcExpt.h"
#endif

#include "hamcRad.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifndef NODICT
ClassImp(hamcRad)
#endif


hamcRad::hamcRad(): did_init(kFALSE)
{
   dE_IntBrehm = 0;
   dE_ExtBrehmIn = 0;
   dE_ExtBrehmOut = 0;
}

hamcRad::~hamcRad()
{
}


#ifdef NOT_STANDALONE
Int_t hamcRad::Init(hamcExpt* expt) {
// Here you want to grab from "expt" the parameters you need
// which depends on experiment and is the same for all events,
// i.e. definition of target

  Float_t tlen   = expt->target->GetLength();
  Float_t trlen  = expt->target->GetRadLength();
  Float_t tgtZ   = expt->target->GetZ();
  Float_t energy = expt->event->beam->GetEnergy();
// The following only works for single-arm.
  Float_t theta  = expt->GetSpectrom(0)->GetScattAngle();

  return Init(energy, theta, tgtZ, trlen, tlen);

}
#endif

Int_t hamcRad::Init(Float_t E, Float_t theta, Float_t z, Float_t rl, Float_t tl) {

// E = beam energy (GeV)
// theta = central angle (degrees)
// Z = density-weighted Z of target
// radlen = fractional radiation length, density-weighted
// len = actual length (meters)

  ldebug = 0;
  Npts = 2000;
  yfact = 1e7;
  nybin = 1000;
  Nslices = 10;  // # slices of target.

  Float_t x1,x2;

  ycell = 0.02*yfact/((Float_t)nybin);
  me=0.511;     // mass electron (MeV)
  alpha=(1./137.);
  pi=3.1415926;

  E0 = E;
  trlen = rl;
  tlen  = tl;

  Float_t theta_rad = theta*3.1415/180;
  qsq = 2*E0*E0*(1-TMath::Cos(theta_rad));
  Float_t msq = (me/1000)*(me/1000);  // GeV^2

  dE = E0/((Float_t)Npts);  // interval

  x1 = exp((-2./3.)*log(z));
  x2 = exp((-1./3.)*log(z));
  Float_t psi = log(1440*x1)/log(183*x2);
  bval = (4./3.)*(1 + (1./9.)*(( (z+1)/(z+psi) ) / (log(183*x2))));

  // This is the equiv. radiator before / after scatt
  tequiv = (alpha/(bval*pi)) * (log(qsq/msq) - 1);

  Setup(0,tequiv);

  Float_t dtgt = trlen/(Float_t(Nslices));

  for (Int_t isl = 0; isl<Nslices; isl++) {

    Float_t tfrac = dtgt * (isl+1);

    Setup(1,tfrac);

  }

  did_init = kTRUE;
  return 1;

}


void hamcRad::Setup(Int_t which, Float_t trl) {

// recall, E0, dE, and bval are global

  Float_t E,Prob,Ptot,Ie;
  Float_t xncell,x1,x2,x3;
  Int_t ncnt,ncell;

  vector<Float_t> radtail;

  E = E0;  // Initialize
  Ptot  = 0; 
  ncnt  = 0;

  for (Int_t i=0; i<Npts; i++) {

    E = E-dE;

    if (E<0) continue;

    x1 = bval*trl/(E0-E);
    x2 = E/E0 + (3./4.)*(((E-E0)/E0)*((E-E0)/E0));
    x3 = exp(bval*trl*log(log(E0/E)));

    if(ldebug==1) cout << "strag. factors "<<i<<"  "<<E<<"  "<<x1<<"  "<<x2<<"  "<<x3<<endl;

    Ie = x1*x2*x3;
    Prob = Ie*dE;
    Ptot += Prob;
    Prob = yfact*Prob;
    xncell = Prob/ycell;
    ncell = ((Int_t)xncell);
    ncnt += ncell;

    if (ncnt > MAXCNT) {
      // this should never happen
      cout << "hamcRad::ERROR:  too many cells !  " <<ncell<<"  "<<ncnt<<endl;
    } else {
      for (Int_t j=0; j<ncell; j++) {
	if (which==0) {
           de_intern.push_back(E);
	} else {
 	   radtail.push_back(E);
	}
      }
    }        
    
  }

  if (which ==0) {
    if(ldebug) cout << "de_intern size "<<de_intern.size()<<endl;
  } else {
    if(ldebug) cout << "de_straggle["<<which<<"] size "<<radtail.size()<<endl;
  }

  if (!which) return;

  de_straggle.push_back(radtail);


}


Int_t hamcRad::LookupIdx(Float_t tl) {

  Float_t dtgt = tlen/(Float_t(Nslices));
  Float_t t1, t2;

  t1=0;

  for (Int_t isl = 0; isl<Nslices; isl++) {

    t2 = dtgt * (isl+1);

    if (tl > t1 && tl <= t2) return isl;

    t1 = t2;

  }

  return Nslices-1;

}


#ifdef NOT_STANDALONE
Int_t hamcRad::Generate(hamcExpt* expt) {
// This routine will be called for each event.
// Main thing you need is the Z location (meters) 
// in target for the main scattering point, to decide
// how much length before / after scattering.

  Float_t ztgt = expt->target->GetZloc();

  return Generate(ztgt);;  

}
#endif

Int_t hamcRad::Generate(Float_t ztgt) {

//  Generate the energy losses in target
//  ztgt = location in target (meters) of scattering
//  (0 = front, and it goes up to tlen)

   dE_IntBrehm = 0;
   dE_ExtBrehmIn = 0;
   dE_ExtBrehmOut = 0;

   Int_t idx, jj;

   if (CheckInit() == -1) return -1;
   
   Float_t x = (de_intern.size()-1)*gRandom->Rndm();
   idx = (Int_t)x;

   dE_IntBrehm = de_intern[idx];

   idx = LookupIdx(ztgt);  // before scattering

   if (idx >= 0 && idx < Nslices) {
      vector<Float_t> radcor = de_straggle[idx];
      x = (radcor.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmIn = radcor[jj];
   }

   idx = LookupIdx(tlen-ztgt);  // after scattering

   if (idx >= 0 && idx < Nslices) {
      vector<Float_t> radcor = de_straggle[idx];
      x = (radcor.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmOut = radcor[jj];
   }

   return 1;

}


Float_t hamcRad::GetDeIntern() {
   
  CheckInit();

  return dE_IntBrehm;

}

Float_t hamcRad::GetDeExternIn() {
   
  CheckInit();

  return dE_ExtBrehmIn;

}

Float_t hamcRad::GetDeExternOut() {
   
  CheckInit();

  return dE_ExtBrehmOut;

}


Int_t hamcRad::CheckInit() {

   if (!did_init) {
    cout <<"hamcRad::ERROR:  Did not initialize the class !"<<endl;
    return -1;
   }

}
