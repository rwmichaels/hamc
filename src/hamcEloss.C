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
#include "hamcSpecHRS.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifndef NODICT
ClassImp(hamcEloss)
#endif

double radtail_tf1(double *x, double *p) { 
// TF1 radtail with x[0] = Eloss
// p[0] = trl (radiation length)
// p[1] = Z (charge)
// p[2] = E0 (initial energy)
  double cut = 1e-10;
  double trl = p[0];
  double z = p[1];
  double E0 = p[2]; 
  double eloss = x[0];
  if (eloss < 0) return 0;
  if (eloss < cut) eloss = cut;
  double E = E0 - eloss;
  double x1,x2,x3,bval,psi;
  x1 = exp((-2./3.)*log(z));
  x2 = exp((-1./3.)*log(z));
  psi = log(1440*x1)/log(183*x2);
  bval = (4./3.)*(1 + (1./9.)*(( (z+1)/(z+psi) ) / (log(183*x2))));
  x1 = bval*trl/eloss;
  x2 = E/E0 + (3./4.)*((eloss/E0)*(eloss/E0));
  x3 = exp(bval*trl*log(log(E0/E)));
  return x1*x2*x3;
}


hamcEloss::hamcEloss(): did_init(kFALSE)
{
   dE_IntBrehm = 0;
   dE_ExtBrehmIn = 0;
   dE_ExtBrehmOut = 0;
   dE_Bsum = 0;
   dE_IonizeIn = 0;
   dE_IonizeOut = 0;

// Possible states, one (and only one) is kTRUE:
//     use_genercone       use_tf1    use_numer
//   (pref. by HAPPEX)   (function)   (numerical)
//
// if use_numer, use_ionize = kFALSE (because its folded in)
// else use_ionize = kTRUE

   use_genercone = kTRUE;  // This is preferred, for now (Jan '09)
   use_ionize = kTRUE;
   use_tf1 = kFALSE;
   use_numer = kFALSE;
}

hamcEloss::~hamcEloss()
{
}


Int_t hamcEloss::Init(hamcExpt* expt) {
// Here you want to grab from "expt" the parameters you need
// which depends on experiment and is the same for all events,
// i.e. definition of target

   hamcStrParser parser;
   
   psi_scale = 0.17;  // default for PREX
   parser.Load(expt->inout->GetStrVect("radcor"));
   if (parser.IsFound("genercone")) {
      cout << "hamcEloss: Using genercone radiative tail "<<endl;
      use_genercone = kTRUE;  
      use_ionize = kTRUE;
      use_tf1 = kFALSE;
      use_numer = kFALSE;
   }   
   if (parser.IsFound("tf1")) {
      cout << "hamcEloss: Using TF1 radiative tail "<<endl;
      use_genercone = kFALSE;  
      use_ionize = kTRUE;
      use_tf1 = kTRUE;
      use_numer = kFALSE;
   }   
   if (parser.IsFound("numer")) {
      cout << "hamcEloss: Using mumerical radiative tail "<<endl;
      use_genercone = kFALSE;  
      use_ionize = kFALSE;  // ionization folded in already
      use_tf1 = kFALSE;
      use_numer = kTRUE;
   }   
   if (parser.IsFound("psi")) {
      psi_scale = parser.GetData(); 
      cout << "hamcEloss:: psi_scale factor = "<<psi_scale<<endl;
   }      

   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dE1", 
	 "dE Brehm in", &dE_ExtBrehmIn, 2000,-0.1,1.2);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "dE2", 
	 "dE Brehm out", &dE_ExtBrehmOut, 2000,-0.1,1.2);
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

// variables belonging to this class
   tlen   = expt->target->GetLength();  // Length (m)
   trlen  = expt->target->GetRadLength();  // Rad len (frac)
   tgtA   = expt->target->GetA();   //  <A> of tgt 
   tgtZ   = expt->target->GetZ();   //  <Z> of tgt 
   tdensity = expt->target->GetDensity();
   E0 = 1;
   if (expt->event->beam) {
      E0 = expt->event->beam->GetE0();
   } else {
      cout << "hamcEloss::WARNING:  No beam, using default energy = "
	  <<E0<<"  GeV"<<endl;
   }
// The following only works for single-arm.
   theta_central = (3.14159/180.0)*
           expt->GetSpectrom(0)->GetScattAngle();


   did_init = kTRUE;

   return InitRad();

}

Int_t hamcEloss::InitRad() {

  if (use_genercone) return OK;  // does not use this

  qsq = 2*E0*E0*(1-TMath::Cos(theta_central));
  Float_t alpha = (1./137.);
  Float_t bval = 4./3.;
  Float_t msq = 2.6112e-7;  // mass electron squared (GeV^2)
// This is the equiv. radiator used for Internal Brehms.
  tequiv = (alpha/(bval*pi)) * (log(qsq/msq) - 1);

  if (use_tf1) {
    Setup_tf1(0, tequiv, tgtZ, E0);  
  }

  Nslices = 10;  // # slices of target.

  Float_t dtgt = trlen/(Float_t(Nslices));

  for (Int_t isl = 1; isl<=Nslices; isl++) {

    Float_t tfrac = dtgt * (isl+1);

    if (use_tf1) {
      Setup_tf1(isl, tfrac, tgtZ, E0);
    } else {
      Setup_numer(isl, tfrac);
    }

  }

  //  Print();

  did_init = kTRUE;
  return 1;

}


void hamcEloss::Setup_numer(Int_t which, Float_t trl) {

// Setup numerical tables that includes internal,
// external, and ionization all in one.
// Recall, E0 is global
// "which" = slices number
// trl = radiation length for this slice size

  Npts = 1000;
  yfact = 10000;  // 1e6 requires MAXCNT=5e6

  Int_t ldebug = 0;

  if (use_genercone) return; // dont need this for genercone version
  if (use_tf1) return;   // error to call this for this state.

  Int_t ncnt,ncell;

  Float_t psi = psi_scale;

  if (psi_scale == 0) {
    psi = 1;
    if (tgtA > 100) psi_scale = 0.17;
  }

  vector<Float_t> radtail;

  Float_t ehi = 1.1*E0;  
  Float_t elo = 0.8*E0;  
  Float_t dE = (ehi-elo)/((Float_t)(Npts));  // GeV
  ncnt = 0;

  for (Int_t i=0; i<Npts; i++) {

    Float_t E = elo + ((Float_t)i)*dE; 
    Float_t eloss = E0-E;
    Float_t lambda = (1000.*eloss/psi) - 0.225;

    Double_t yval = 0;

    if (lambda > 0 && lambda < 800) {
       yval = GenLandHi(lambda, (trl+0.5*tequiv));
    } 
    if (lambda > -10 && lambda <= 0) {
       yval = GenLandLo(lambda, (trl+0.5*tequiv));
    }

    yval = yval*yfact; 

    ncell = ((Int_t)yval);
    ncnt += ncell;

    if (ldebug && lambda > -6 && lambda < 500) cout << "tail "<<i<<"  "<<trl<<"  "<<E<<"  "<<eloss<<"  "<<lambda<<"  "<<yval<<"  "<<ncell<<"  "<<ncnt<<endl;

    if (ncnt > MAXCNT) {
      // this should never happen
      cout << "hamcEloss::ERROR:  too many cells !  " <<ncell<<"  "<<ncnt<<endl;
    } else {

      for (Int_t j=0; j<ncell; j++) {
 	   radtail.push_back(eloss);
      }

    }

  }

  Enumer.push_back(radtail);

}

void hamcEloss::Setup_tf1(Int_t which, Double_t rad, Double_t Ztgt, Double_t Ene) {
// Setup the radiative tail based on a TF1 radtail_tf1.
// "which" flag: 0=internal, >0 = slice#.
// radiation length = rad
// charge of target material = Ztgt
// initial electron energy = Ene

   cout << "TF1 Set up "<<which<<"  "<<rad<<"  "<<Ztgt<<"  "<<Ene<<endl;
   char ctf1name[80];
   sprintf(ctf1name,"distrFunc%d",which);
   fdistr.push_back(new TF1(ctf1name,radtail_tf1,0.,Ene,3));
   Double_t par[3];
   par[0] = rad; 
   par[1] = Ztgt;
   par[2] = Ene;
   Int_t idx = fdistr.size()-1;
   fdistr[idx]->SetParameters(par); 
   cout << "TF1  fdistr["<<idx<<"] =  "<<fdistr[idx]<<"   rad "<<rad<<endl;

}

void hamcEloss::Print() {

  cout << "hamcEloss Print: "<<endl;
  cout << "Enumer size "<<Enumer.size()<<endl;
  if ((Int_t)Enumer.size() != Nslices) {
    cout << "Nslices = "<<Nslices<<"  inconsistent "<<endl;
    return;
  }
  
  vector<Float_t> radtail;

  for (Int_t jj = 0; jj<Nslices; jj++) {
    cout << endl << "Brehms. slice # "<<jj<<endl;
    radtail = Enumer[jj];
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

    if (tl > t1 && tl <= t2) return isl+1;

    t1 = t2;

  }

  return 0;

}


Int_t hamcEloss::Generate(hamcExpt* expt) {
// This routine will be called for each event.
// Main thing you need is the Z location (meters) 
// in target for the main scattering point, to decide
// how much length before / after scattering.

  dE_IntBrehm = 0;  // initialize
  dE_ExtBrehmIn = 0;
  dE_ExtBrehmOut = 0;
  dE_Bsum = 0;
  dE_IonizeIn = 0;
  dE_IonizeOut = 0;
  radin = expt->target->GetRadIn();   // variable of the class
  radout = expt->target->GetRadOut(); //  "    "     "
  
  if (use_genercone) {
     Generate_gcone();
     return 1;
  } 

  if (use_tf1) {
     Generate_tf1();
     return 1;
  }

  if (use_numer) {
     GenerateNumer();
     return 1;
  }

  if (!use_numer) GenerateDeDx(expt);  

  return OK;

}


Float_t hamcEloss::GetDeDx(Float_t radlen, Int_t where) {
 
  if (where != ICOLLIM2) return 0;  // only for 2nd collimator -- Be plug

  // For Be the fraction is 0.1023
  // For C12 it is 0.0733
  Float_t dE = 0.0733*radlen;  // GeV per fraction
 // e.g. 1.5 cm Be = 4.3 % RL = 0.043 -> 0.0044 GeV = 4.4 MeV

  return dE;

}


Int_t hamcEloss::GenerateDeDx(hamcExpt *expt) {
// Generate the de/dx energy losses due to ionization
// in the target.  For degrader collimator like Be plug
// use GetDeDx(Float_t radlen, Int_t where);
// Model here uses no straggling.
// Note, this is already built into the numerical results (if 'use_numer')

  dE_IonizeIn = 0; 
  dE_IonizeOut = 0;

  if ( !use_ionize ) return 0;
  if ( use_numer ) return 0;   // already built in

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

  return OK;

}

Int_t hamcEloss::Generate_gcone() {

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

Int_t hamcEloss::Generate_tf1() {

// Version uses a TF1 random number generator.
// Input: zpos = position (meters) in target.
// goes from -tlen/2 to +tlen/2.

   if (fdistr.size() == 0) return -1;

   Int_t idx;

   dE_IntBrehm = fdistr[0]->GetRandom();

   idx = LookupIdx(radin);  // before scattering

   if (idx > 0 && idx < (Int_t)fdistr.size()) {
        dE_ExtBrehmIn = fdistr[idx]->GetRandom();
   }

   idx = LookupIdx(radout);  // after scattering

   if (idx > 0 && idx < (Int_t)fdistr.size()) {
        dE_ExtBrehmOut = fdistr[idx]->GetRandom();
   }

   dE_Bsum = dE_IntBrehm + dE_ExtBrehmIn + dE_ExtBrehmOut;

   return 1;
}


Int_t hamcEloss::GenerateNumer() {

// This version uses numerical radiatiave tails 
// generated in Setup.
// Input: zpos = position (meters) in target.
// goes from -tlen/2 to +tlen/2.

   Int_t idx, jj;
   Float_t x;

   if (CheckInit() == -1) return -1;
   
   idx = LookupIdx(radin);  // before scattering

   vector<Float_t> radtail;  

   if (idx >= 0 && idx < Nslices) {
      radtail = Enumer[idx];
      if (radtail.size() == 0) {
        cout << "hamcEloss: WARNING: no radtail defined for slice "<<idx<<endl;
        return 0;
      }
      x = (radtail.size()-1.)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmIn = radtail[jj];
   }

   idx = LookupIdx(radout);  // after scattering

   if (idx >= 0 && idx < Nslices) {
      radtail = Enumer[idx];
      if (radtail.size() == 0) {
        cout << "hamcEloss: WARNING: no radtail defined for slice "<<idx<<endl;
        return 0;
      }
      x = (radtail.size()-1)*gRandom->Rndm();
      jj = (Int_t)x; 
      dE_ExtBrehmOut = radtail[jj];
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

  // Model for stopping power from PDG fit.
  // distance = distance through target (meters)
  // No straggling (Landau tail) in this version.
  // Ideally, need better treatment of composite target.

  if (Znuc <= 0) {
    cout << "hamcEloss::ERROR: trying to calc. dE/dx for neutral ?"<<endl;
    return 0;
  }

  Float_t dedx = 2.35 - 0.28*TMath::Log(Znuc);  // MeV g^-1 cm^2

  // exception for hydrogen
  if (Znuc == 1) dedx = 3.97;  // MeV g^-1 cm^2

  return dedx;

}

   
Double_t hamcEloss::GenLandLo(Float_t lambda, Float_t trl) {
   
  Double_t rmax = 1;
  Double_t rmin = 0.001;
  Double_t xpts = 5000;
  Int_t npts = (Int_t)xpts;
  Double_t dr = (rmax - rmin) / xpts;

  Double_t dd = TMath::Exp(-1.0*(lambda + 1));

  Double_t xfact1 = (1./pi)*TMath::Exp(-1.0*dd*(1 + (trl/dd)*TMath::Log(dd)));

  Double_t result = 0;

  for (Int_t i = 0; i < npts; i++) {

    Double_t r = rmin + (rmax-rmin)*((Double_t)i)/xpts;

    result += dr * xfact1 * (
	 TMath::Exp (
          0.5*(dd-trl) * TMath::Log(1.0 + (r/dd)*(r/dd)) 
	  - dd*TMath::ATan(r/dd) )
	 * TMath::Cos (r*(0.5*TMath::Log(1.0 + (r/dd)*(r/dd)) -1.0))
         + (dd-trl)*TMath::ATan(r/dd));

  }

  return result;

}


Double_t hamcEloss::GenLandHi(Float_t lambda, Float_t trl) {

  Double_t rmax = 1;
  Double_t rmin = 0.001;
  Double_t xpts = 5000;
  Int_t npts = (Int_t)xpts;
  Double_t dr = (rmax - rmin) / xpts;

  Double_t result = 0;

  for (Int_t i = 0; i < npts; i++) {

    Double_t r = rmin + (rmax-rmin)*((Double_t)i)/xpts;

    result += (dr/pi) * (TMath::Sin(pi*(trl + r)) * 
	       TMath::Exp(-1.0*lambda*r - (trl+r)*TMath::Log(r)));

  }

  return result;
  
}
