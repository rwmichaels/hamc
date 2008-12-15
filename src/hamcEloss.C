//  hamcEloss   -- Energy losses 
//  Internal and external Brehmstrahlung
//  Ionization losses.
//  Member of hamcPhysics

//  R. Michaels  Nov 2008


#include "hamcExpt.h"
#include "hamcTarget.h"
#include "hamcEvent.h"
#include "hamcBeam.h"
#include "hamcEloss.h"
#include "hamcInout.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifndef NODICT
ClassImp(hamcEloss)
#endif


hamcEloss::hamcEloss(): did_init(kFALSE)
{
   dE_IntBrehm = 0;
   dE_ExtBrehmIn = 0;
   dE_ExtBrehmOut = 0;
   dE_Bsum = 0;
   dE_IonizeIn = 0;
   dE_IonizeOut = 0;
   use_genercone = kTRUE;  // For tests. (might become permanent)
   use_ionize = kTRUE;
}

hamcEloss::~hamcEloss()
{
}


Int_t hamcEloss::Init(hamcExpt* expt) {
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
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dIonin", 
	 "dE Ionization loss in", &dE_IonizeIn, 1000,-0.02,0.1);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dIonout", 
	 "dE Ionization loss out", &dE_IonizeOut, 1000,-0.02,0.1);

   expt->inout->AddToNtuple("dE1",&dE_ExtBrehmIn);
   expt->inout->AddToNtuple("dE2",&dE_ExtBrehmOut);
   expt->inout->AddToNtuple("dEi",&dE_IntBrehm);
   expt->inout->AddToNtuple("dIin",&dE_IonizeIn);
   expt->inout->AddToNtuple("dIout",&dE_IonizeOut);


   tlen   = expt->target->GetLength();
   trlen  = expt->target->GetRadLength();
   tgtZ   = expt->target->GetZ();
   tdensity = expt->target->GetDensity();
   Float_t energy = 1;
   if (expt->event->beam) {
      energy = expt->event->beam->GetE0();
   } else {
      cout << "hamcEloss::WARNING:  No beam, using default energy = "
	  <<energy<<"  GeV"<<endl;
   }
// The following only works for single-arm.
   Float_t theta  = expt->GetSpectrom(0)->GetScattAngle();

   return InitRad(energy, theta, tgtZ, trlen, tlen);

}

Int_t hamcEloss::InitRad(Float_t E, Float_t theta, Float_t z, Float_t rl, Float_t tl) {

// E = beam energy (GeV)
// theta = central angle (degrees)
// Z = density-weighted Z of target
// radlen = fractional radiation length, density-weighted
// len = actual length (meters)

  Npts = 50000;
  yfact = 2e5;
  Nslices = 10;  // # slices of target.

  if (ldebug) cout << "hamcEloss init "<<E<<"  "<<theta<<"  "<<z<<"  "<<rl<<"  "<<tl<<endl;  

  Float_t x1,x2;

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

  //  Print();

  did_init = kTRUE;
  return 1;

}


void hamcEloss::Setup(Int_t which, Float_t trl) {

// recall, E0, dE, and bval are global

  Float_t E,eresol,Prob,XProb,Ptot,Ie;
  Float_t x1,x2,x3,phi,plo,phltot,pcut;
  Int_t ncnt,ncell;

  if (use_genercone) return; // dont need this for genecone version

  vector<Float_t> radtail;

  E = E0;  // Initialize
  eresol = fracresol*E0;
  Ptot  = 0; 
  ncnt  = 0;

  for (Int_t i=0; i<Npts; i++) {

    E = E-dE;

    if (E<0) continue;

    x1 = bval*trl/(E0-E);
    x2 = E/E0 + (3./4.)*(((E-E0)/E0)*((E-E0)/E0));
    x3 = exp(bval*trl*log(log(E0/E)));

    if(ldebug==1) cout << "rad. factors "<<i<<"  "<<E<<"  "<<x1<<"  "<<x2<<"  "<<x3<<endl;

    Ie = x1*x2*x3;
    Prob = Ie*dE;
    if ((E0-E) > eresol) {
      phi += Prob;
    } else {
      plo += Prob;
    }
    Ptot += Prob;
    XProb = yfact*Prob;
    ncell = ((Int_t)XProb);
    ncnt += ncell;

    if (ncnt > MAXCNT) {
      // this should never happen
      cout << "hamcEloss::ERROR:  too many cells !  " <<ncell<<"  "<<ncnt<<endl;
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

  phltot = phi+plo;  // This should be ~1, but may fall a bit short.
  pcut = 0;
  if (phltot != 0) pcut = plo/phltot;   

  if (which ==0) {
    cut_intern = pcut;
    if(ldebug) cout << "Eintern size "<<Eintern.size()<<endl;
  } else {
    if(ldebug) cout << "Eextern["<<which<<"] size "<<radtail.size()<<endl;
  }

  if (!which) return;

  Eextern.push_back(radtail);
  cut_extern.push_back(pcut);

}

void hamcEloss::Print() {

  cout << "hamcEloss Print: "<<endl;
  cout << "Eintern size "<<Eintern.size()<<endl;
  cout << "Straggle # slices "<<Eextern.size()<<endl;
  cout << "Ecut for int. Brehm "<<cut_intern<<endl;
  cout << "Ecut for ext. Brehm by slice ";
  for (Int_t i = 0; i<(Int_t)cut_extern.size(); i++) {
    cout << "  cut["<<i<<"] = "<<cut_extern[i];
  }
  cout << endl;
  
  cout << endl << "Internal Brehms. tail"<<endl;

  for (Int_t i = 0; i < (Int_t)Eintern.size(); i++) {
    cout << "  E("<<i<<") = "<<Eintern[i];
    if (i > 0 && (i%8==0)) cout << endl;
  }

  vector<Float_t> radtail;

  for (Int_t jj = 0; jj<Nslices; jj++) {
    cout << endl << "External Brehms. # "<<jj<<endl;
    radtail = Eextern[jj];
    for (Int_t i = 0; i < (Int_t)radtail.size(); i++) {
      cout << "  E("<<i<<")  = "<<radtail[i];
      if (i > 0 && (i%8==0)) cout << endl;
    }
  }


}


Int_t hamcEloss::LookupIdx(Float_t tl) {

  Float_t dtgt = trlen/(Float_t(Nslices));
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


Int_t hamcEloss::Generate(hamcExpt* expt) {
// This routine will be called for each event.
// Main thing you need is the Z location (meters) 
// in target for the main scattering point, to decide
// how much length before / after scattering.

  Float_t radin = expt->target->GetRadIn();
  Float_t radout = expt->target->GetRadOut();
  
  if (use_genercone) {
    GenerateRad(radin, radout);
  } else {
    GenerateRad(expt->target->GetZScatt());
  }

  GenerateDeDx(expt);

  return OK;

}

Int_t hamcEloss::GenerateDeDx(hamcExpt *expt) {
// Generate the de/dx energy losses due to ionization.
// Model uses no straggling.

  dE_IonizeIn = 0; 
  dE_IonizeOut = 0;

  if ( !use_ionize ) return 0;

  Float_t zmtl, density, zlen;

  Int_t num_mtl = expt->target->GetNumMtl();
  Int_t mtl_idx = expt->target->GetMtlIndex();

  if (mtl_idx < 0 || mtl_idx >= num_mtl) {
    cout << "hamcEloss::ERROR: unexpected material index"<<endl;
    return 0;
  }

// Materials prior to scatter point
  for (Int_t islab = 0; islab<mtl_idx; islab++) {

    density = expt->target->GetMtlDensity(islab);
    zlen = expt->target->GetMtlLen(islab);    
    zmtl = expt->target->GetMtlZ(islab);

  // density is g/cm^2 and dist in meters, so must convert
  // 100 cm/meter * dist(meters) * 
  // 0.001 GeV/MeV  -> overall factor of 0.1

    dE_IonizeIn += 0.1*dedx_eloss(zmtl)*density*zlen;
   
  }

// Material from which we scattered:

    density = expt->target->GetMtlDensity(mtl_idx);
    zlen = expt->target->GetLenInMtl();
    zmtl = expt->target->GetMtlZ(mtl_idx);
    dE_IonizeIn += 0.1*dedx_eloss(zmtl)*density*zlen;

// Now going out from scattering point
    zlen = expt->target->GetLenInMtl();
    dE_IonizeOut += 0.1*dedx_eloss(zmtl)*density*zlen;

// Materials after the scatter point
  for (Int_t islab = num_mtl-1; islab>mtl_idx; islab--) {

    density = expt->target->GetMtlDensity(islab);
    zlen = expt->target->GetMtlLen(islab);    
    zmtl = expt->target->GetMtlZ(islab);
    dE_IonizeOut += 0.1*dedx_eloss(zmtl)*density*zlen;
   
  }

}



Int_t hamcEloss::GenerateRad(Float_t radin, Float_t radout) {

//  This version uses gener_cone methods
//  Generate the radiative losses in target
//  radin = # radiation lengths coming in to scatter point
//  radout = # rad  "    "  out "  "

   dE_IntBrehm = gener_radlossint(E0,0.5*tequiv);
   dE_ExtBrehmIn = gener_radlossext(E0, radin);
   dE_ExtBrehmOut = gener_radlossext(E0, radout);
   dE_Bsum = dE_IntBrehm + dE_ExtBrehmIn + dE_ExtBrehmOut;

     //     cout << "Using Gener cone "<<E0<<"  "<<tequiv<<"  "<<"  "<<"  "<<radin<<"  "<<radout<<" / "<<endl<<dE_IntBrehm<<"  "<<dE_ExtBrehmIn<<"  "<<dE_ExtBrehmOut<<"  "<<dE_Bsum<<endl;

     return 1;

}

Int_t hamcEloss::GenerateRad(Float_t zpos) {

// This version uses radiatiave tails generated in Setup.
// Input: zpos = position (meters) in target.
// goes from -tlen/2 to +tlen/2.

   Int_t idx, jj;
   Float_t prob;

   if (CheckInit() == -1) return -1;
   
   Float_t x = (Eintern.size()-1)*gRandom->Rndm();
   idx = (Int_t)x;
   if (idx < 0 || idx > (Int_t)Eintern.size()) return -1;


   prob = gRandom->Rndm();
   if (prob > cut_intern) {  // if there was an Brehm. event
     dE_IntBrehm = E0 - Eintern[idx];
   }

   Float_t xin = zpos + tlen/2;
   Float_t rin = xin*trlen/tlen;

   idx = LookupIdx(rin);  // before scattering

   //   cout << "Stuff1 "<<tlen<<"  "<<trlen<<"  "<<Nslices<<"  "<<xin<<"   "<<rin<<"  "<<idx<<endl;

   vector<Float_t> radtail;  

   if (idx >= 0 && idx < Nslices) {
      radtail = Eextern[idx];
      x = (radtail.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      prob = gRandom->Rndm();
      if (prob > cut_extern[idx]) {// if there was an Brehm. event
        dE_ExtBrehmIn = E0 - radtail[jj];
      }
   }

   Float_t xout = tlen/2 - zpos;
   Float_t rout = xout*trlen/tlen;

   idx = LookupIdx(rout);  // after scattering

   //   cout << "Stuff2 "<<rout<<"  "<<idx<<endl;

   if (idx >= 0 && idx < Nslices) {
      radtail = Eextern[idx];
      x = (radtail.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmOut = E0 - radtail[jj];
   }


   dE_Bsum = dE_IntBrehm + dE_ExtBrehmIn + dE_ExtBrehmOut;

   return 1;


}


Float_t hamcEloss::GetDeIntern() {
   
  CheckInit();

  return dE_IntBrehm;

}

Float_t hamcEloss::GetDeExternIn() {
   
  CheckInit();

  return dE_ExtBrehmIn;

}

Float_t hamcEloss::GetDeExternOut() {
   
  CheckInit();

  return dE_ExtBrehmOut;

}


Int_t hamcEloss::CheckInit() {

   if (!did_init) {
    cout <<"hamcEloss::ERROR:  Did not initialize the class !"<<endl;
    return -1;
   }

   return 1;
}


Float_t hamcEloss::gener_radlossint(Float_t k, Float_t nu) {

/* Taken from gener_cone MC, modified by R. Michaels */  

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


Float_t hamcEloss::gener_radlossext(Float_t k, Float_t fracrl) {

/* Taken from gener_cone MC, modified by R.M. */
/* For now we ignore the subtleties of target composition, treat it like
   one big mix of an effective target length; will fix this later. - Bob */

/* L.Van Hoorebeke, University of Gent, e-mail: Luc.VanHoorebeke@UGent.be, luc@inwfsun1.rug.ac.be
  This function generates external bremsstrahlung in the target
  the distribution used contains multiple emission effects
  k: electron momentum
  fracrl = fraction of radiation lengths distance, travelled through material.

*/

  Float_t cut,Ekin,bt,prob,prob_sample;
  Float_t sample,xtry,env,value,ref;

/* Initialisation of lower limit of bremsstrahlung (1 keV) */
  cut=0.000001;

  Ekin = k-Me;

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


Float_t hamcEloss::dedx_eloss(Float_t Znuc) {

  // Model for stopping power from PDG fit
  // distance = distance through target (meters)
  // No straggling in this version.

  // FIXME: need better treatment of composite target.

  Float_t dedx = 2.35 - 0.28*TMath::Log(Znuc);  // MeV g^-1 cm^2

  // exception for hydrogen
  if (tgtZ == 1) dedx = 3.97;  // MeV g^-1 cm^2

  return dedx;

}

   
