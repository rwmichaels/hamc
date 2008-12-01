//  hamcRad   -- Internal and external Brehmstrahlung
//  Member of hamcPhysics

//  R. Michaels  Nov 2008


#include "hamcExpt.h"
#include "hamcTarget.h"
#include "hamcEvent.h"
#include "hamcBeam.h"
#include "hamcRad.h"
#include "hamcInout.h"
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
   dE_Bsum = 0;
   use_genercone = kTRUE;  // For tests
}

hamcRad::~hamcRad()
{
}


Int_t hamcRad::Init(hamcExpt* expt) {
// Here you want to grab from "expt" the parameters you need
// which depends on experiment and is the same for all events,
// i.e. definition of target


   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dE1", 
	 "dE Brehm in", &dE_ExtBrehmIn, 1000,-0.1,1.2);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dE2", 
	 "dE Brehm out", &dE_ExtBrehmOut, 1000,-0.1,1.2);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dEi", 
	 "dE Brehm intern", &dE_IntBrehm, 1000,-0.1,1.2);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dEsum", 
	 "dE Sum Brehm", &dE_Bsum, 1000,-0.1,1.2);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dE2i", 
	 "dE Brehm correl (2,i)", 
              &dE_IntBrehm, 100,-0.05,0.5,
              &dE_ExtBrehmOut, 100,-0.05,0.5);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dE12", 
	 "dE Brehm correl (1,2)", 
              &dE_ExtBrehmIn, 100,-0.05,0.5,
              &dE_ExtBrehmOut, 100,-0.05,0.5);

   expt->inout->AddToNtuple("dE1",&dE_ExtBrehmIn);
   expt->inout->AddToNtuple("dE2",&dE_ExtBrehmOut);
   expt->inout->AddToNtuple("dEi",&dE_IntBrehm);



  Float_t tlen   = expt->target->GetLength();
  Float_t trlen  = expt->target->GetRadLength();
  Float_t tgtZ   = expt->target->GetZ();
  Float_t energy = 1;
  if (expt->event->beam) {
      energy = expt->event->beam->GetE0();
  } else {
    cout << "hamcRad::WARNING:  No beam, using default energy = "
	 <<energy<<"  GeV"<<endl;
  }
// The following only works for single-arm.
  Float_t theta  = expt->GetSpectrom(0)->GetScattAngle();



  return Init(energy, theta, tgtZ, trlen, tlen);

}

Int_t hamcRad::Init(Float_t E, Float_t theta, Float_t z, Float_t rl, Float_t tl) {

// E = beam energy (GeV)
// theta = central angle (degrees)
// Z = density-weighted Z of target
// radlen = fractional radiation length, density-weighted
// len = actual length (meters)

  Npts = 2000;
  yfact = 1e7;
  nybin = 1000;
  Nslices = 10;  // # slices of target.

  if (ldebug) cout << "hamcRad init "<<E<<"  "<<theta<<"  "<<z<<"  "<<rl<<"  "<<tl<<endl;  

  Float_t x1,x2;

  ycell = 0.01*yfact/((Float_t)nybin);
  me=0.511;     // mass electron (MeV)
  alpha=(1./137.);
  pi=3.1415926;

  E0 = E;
  trlen = rl;
  tlen  = tl;

  Float_t theta_rad = theta*pi/180;
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

  if (use_genercone) return; // dont need this for genecone version

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
           Eintern.push_back(E);
	} else {
 	   radtail.push_back(E);
	}
      }
    }        
    
  }

  if (which ==0) {
    if(ldebug) cout << "Eintern size "<<Eintern.size()<<endl;
  } else {
    if(ldebug) cout << "Estraggle["<<which<<"] size "<<radtail.size()<<endl;
  }

  if (!which) return;

  Estraggle.push_back(radtail);


}

void hamcRad::Print() {

  cout << "hamcRad Print: "<<endl;
  cout << "Eintern size "<<Eintern.size()<<endl;
  cout << "Straggle # slices "<<Estraggle.size()<<endl;
  
  cout << endl << "First part of internal "<<endl;

  for (Int_t i = 0; i < (Int_t)Eintern.size(); i++) {
    cout << "  E("<<i<<") = "<<Eintern[i];
    if (i > 0 && (i%8==0)) cout << endl;
  }

  vector<Float_t> radtail;

  for (Int_t jj = 0; jj<Nslices; jj++) {
    cout << endl << "First part of external # "<<jj<<endl;
    radtail = Estraggle[jj];
    for (Int_t i = 0; i < (Int_t)radtail.size(); i++) {
      cout << "  E("<<i<<")  = "<<radtail[i];
      if (i > 0 && (i%8==0)) cout << endl;
    }
  }


}


Int_t hamcRad::LookupIdx(Float_t tl) {

  Float_t dtgt = tlen/(Float_t(Nslices));
  Float_t t1, t2;

  t1=0;
  if (tl < 0) return 0;

  for (Int_t isl = 0; isl<Nslices; isl++) {

    t2 = dtgt * (isl+1);

    if (tl > t1 && tl <= t2) return isl;

    t1 = t2;

  }

  return 0;

}


Int_t hamcRad::Generate(hamcExpt* expt) {
// This routine will be called for each event.
// Main thing you need is the Z location (meters) 
// in target for the main scattering point, to decide
// how much length before / after scattering.

  Float_t ztgt = expt->target->GetZScatt();

  return Generate(ztgt);

}

Int_t hamcRad::Generate(Float_t ztgt) {

//  Generate the energy losses in target
//  ztgt = location in target (meters) of scattering
//  (0 = front, and it goes up to tlen)
// Add 1/2 tgt len since generated ztgt normally 
// starts at -tlen/2 and goes to +tlen/2.

   ztgt = ztgt + tlen/2;

   dE_IntBrehm = 0;
   dE_ExtBrehmIn = 0;
   dE_ExtBrehmOut = 0;

   if (use_genercone) {  // Test using gener_cone version.

     dE_IntBrehm = gener_radlossint(E0,0.5*tequiv);
     dE_ExtBrehmIn = gener_radlossext(E0, ztgt);
     dE_ExtBrehmOut = gener_radlossext(E0, tlen-ztgt);
     dE_Bsum = dE_IntBrehm + dE_ExtBrehmIn + dE_ExtBrehmOut;

     //     cout << "Using Gener cone "<<E0<<"  "<<tequiv<<"  "<<trlen<<"  "<<ztgt<<"  "<<tlen<<" / "<<endl<<dE_IntBrehm<<"  "<<dE_ExtBrehmIn<<"  "<<dE_ExtBrehmOut<<"  "<<dE_Bsum<<endl;


     return 1;

   }


   Int_t idx, jj;

   if (CheckInit() == -1) return -1;
   
   Float_t x = (Eintern.size()-1)*gRandom->Rndm();
   idx = (Int_t)x;

   if (idx < 0 || idx > Eintern.size()) return -1;

   dE_IntBrehm = E0 - Eintern[idx];

   idx = LookupIdx(ztgt);  // before scattering

   vector<Float_t> radtail;  

   if (idx >= 0 && idx < Nslices) {
      radtail = Estraggle[idx];
        x = (radtail.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmIn = E0 - radtail[jj];
   }

   //   cout << "\nidx before  "<<ztgt<<"  "<<idx<<"  "<<dE_ExtBrehmIn<<"  "<<radtail.size()<<endl;

   idx = LookupIdx(tlen-ztgt);  // after scattering

   if (idx >= 0 && idx < Nslices) {
      radtail = Estraggle[idx];
      x = (radtail.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmOut = E0 - radtail[jj];
   }

   //   cout << "idx after  "<<tlen-ztgt<<"  "<<idx<<"  "<<dE_ExtBrehmOut<<"  "<<radtail.size()<<endl;

   dE_Bsum = dE_IntBrehm + dE_ExtBrehmIn + dE_ExtBrehmOut;

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

   return 1;
}


Float_t hamcRad::gener_radlossint(Float_t k, Float_t nu) {

/* Taken from gener_cone MC, modified by R.M. */  

/* L.Van Hoorebeke, University of Gent, e-mail: Luc.VanHoorebeke@UGent.be
   This function generates internal radiation the distribution used
   contains multiple emission effects
   k: electron momentum (GeV)
   nu: equivalent radiator length in units radiation length
       
       this is actually half the equivalent radiator length
       because 1/2 placed before and 1/2 placed after proton
*/

  Float_t cut,Ekin,prob,prob_sample,sample;

/* Initialisation of lower limit of bremsstrahlung (1 keV) */
  cut = 0.000001;

  Ekin = k-Me;

/* Calculation of probability to have internal radiation effect above 1 keV. */
  prob = 1.-pow(cut/Ekin,nu);

  prob_sample = (Float_t)gRandom->Rndm();  /* Random sampling */

  if (prob_sample > prob) return 0.;

/* bremsstrahlung has taken place! Generate photon energy */
  sample = gRandom->Rndm();
  return Ekin*pow(sample*prob + pow(cut/Ekin,nu),1./nu);

}


Float_t hamcRad::gener_radlossext(Float_t k, Float_t dist_tgt) {

/* Taken from gener_cone MC, modified by R.M. */
/* For now we ignore the subtleties of target composition, treat it like
   one big mix of an effective target length; will fix this later. - Bob */

/* L.Van Hoorebeke, University of Gent, e-mail: Luc.VanHoorebeke@UGent.be, luc@inwfsun1.rug.ac.be
  This function generates external bremsstrahlung in the target
  the distribution used contains multiple emission effects
  k: electron momentum
  dist_tgt:distance travelled through material with rad length rd_tgt (cm)

*/

  Float_t cut,Ekin,fracrl,bt,prob,prob_sample;
  Float_t sample,xtry,env,value,ref;

/* Initialisation of lower limit of bremsstrahlung (1 keV) */
  cut=0.000001;

  Ekin = k-Me;
/* Determination of total radiation lenght fraction travelled through */
  fracrl = trlen*dist_tgt/tlen;

  bt=fracrl*4./3.;

/* Calculation of probability to have bremsstrahlung effect above 1 keV */
  prob = 1.- pow(cut/Ekin,bt) - bt/(bt+1.)*(1.- pow(cut/Ekin,bt+1.))
        + 0.75*bt/(2.+bt)*(1.- pow(cut/Ekin,bt+2.));
  prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+PI*PI/6.)); /* Gamma function */

  prob_sample = gRandom->Rndm();
  if (prob_sample > prob) return 0.;

/* Bremsstrahlung has taken place! Generate photon energy with sample and reject,
   using 1/x as envelope */

  do {
   sample = gRandom->Rndm();
   xtry = cut*pow(Ekin/cut,sample);

   env = 1./xtry;
   value = 1./xtry*(1.-xtry/Ekin+0.75*pow(xtry/Ekin,2))*pow(xtry/Ekin,bt);

   sample = gRandom->Rndm();
   ref = value/env;
  } while (sample > ref);

  return xtry;

}
