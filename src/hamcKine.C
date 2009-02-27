//  hamcKine   -- Kinematics generator.
//  Member of hamcPhysics
//  Kinematics are generated uniform in phase space.
//  Works for elastic and DIS.

//  R. Michaels  Nov 2008


#include "hamcKine.h"
#include "hamcExpt.h"
#include "hamcPhysics.h"
#include "hamcEloss.h"
#include "hamcEvent.h"
#include "hamcTarget.h"
#include "hamcInout.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifndef NODICT
ClassImp(hamcKine)
#endif


hamcKine::hamcKine(): did_init(kFALSE)
{
  xbjlo = 0.01;  // Default cuts for DIS
  xbjhi = 0.99;
  qsqlo = 1;
  wsqlo = 4;
  Clear();
}

hamcKine::~hamcKine()
{
}

void hamcKine::Clear() {

  energy = 0; theta = 0;  phi = 0;  qsq = 0;  
  eprime = 0;  y = 0;  x = 0;  

}

Int_t hamcKine::Init(hamcExpt* expt) {
// Here you want to grab from "expt" the parameters you need
// which depends on experiment and is the same for all events,
// i.e. definition of target

  Float_t thetamin, thetamax, phimin, phimax, emin, emax;

  string process = expt->physics->GetProcess();
  Float_t ebeam = expt->event->beam->GetEnergy();
// theta_central, converted to radians, only for single-arm.
  Float_t theta = PI * expt->GetSpectrom(0)->GetScattAngle() / 180;
  Float_t masst = expt->target->GetMass();

// Ranges of angles for the simulation

// Defaults if not defined in input file
   
   thetamin = 0.1*theta;
   thetamax = 3.0*theta;
   phimin = 0;
   phimax = PI;

// Obtain the ranges from the input file if it exists
   
   hamcStrParser parser;
   parser.Load(expt->inout->GetStrVect("out_angles"));
   if (parser.IsFound("thetamin")) {
    parser.Print();
    thetamin = PI*parser.GetData()/180; 
   }    
   if (parser.IsFound("thetamax")) {
     thetamax = PI*parser.GetData()/180; 
   }
   if (parser.IsFound("phimin")) {
     phimin = PI*parser.GetData()/180; 
   }
   if (parser.IsFound("phimax")) {
     phimax = PI*parser.GetData()/180; 
   }

//   dP0_iter = 0;
//   dtheta_iter = 0;
//   dphi_iter = 0;

//   parser.Load(expt->inout->GetStrVect("kick:track"));
//   parser.Print();
//   if (parser.IsFound("P0")) {
//     if (expt->inout->numiter > 1) dP0_iter = parser.GetData();  
//     cout << "Will iterate P0 by "<<100*dP0_iter<<" %"<<endl;
//   }   
//   if (parser.IsFound("theta")){
//     if (expt->inout->numiter > 1) dtheta_iter = parser.GetData();
//     cout << "Will iterate theta by "<<100*dtheta_iter<<" % of its range"<<endl;
//   }
//   if (parser.IsFound("phi")){
//     if (expt->inout->numiter > 1) dphi_iter = parser.GetData();
//     cout << "Will iterate phi by "<<100*phi<<" % of its range"<<endl;
//   }


// For DIS, obtain the min and max output energies
// and the cuts on x, qsq, and wsq that define DIS.
   emin = 0;  emax = 0;
   parser.Load(expt->inout->GetStrVect("dis_setup"));
   if (parser.IsFound("emin")) {
     emin = parser.GetData();
   }
   if (parser.IsFound("emax")) {
     emax = parser.GetData();
   }
   if (parser.IsFound("xbjlo")) {
     xbjlo = parser.GetData();
   }
   if (parser.IsFound("xbjhi")) {
     xbjhi = parser.GetData();
   }
   if (parser.IsFound("qsqlo")) {
     qsqlo = parser.GetData();
   }
   if (parser.IsFound("wsqlo")) {
     wsqlo = parser.GetData();
   }

   cout << "hamcKine:: angles ranges: "<<thetamin<<"  "<<thetamax<<"  "<<phimin<<"  "<<phimax<<endl;
   if (process == "dis" || process == "DIS") {
     cout << "eprime range   "<<emin<<"  "<<emax<<endl;
     cout << "DIS cuts   "<<xbjlo<<"  "<<xbjhi<<"   "<<qsqlo<<"   "<<wsqlo<<endl;
   }

     
   return Init(process, ebeam, theta, masst, thetamin, thetamax, phimin, phimax, emin, emax);

}

Int_t hamcKine::Init(string proc, Float_t eb, Float_t theta, 
		Float_t masst,
		Float_t thmi, Float_t thma, Float_t phmi, Float_t phma,
		Float_t epmi, Float_t epma) {
// energy and mass in GeV
// angles in radians


    cout << "hamcKine: generated process = ";
    iproc = proc_undef;
    if (proc == "elastic" || proc == "Elastic") {
        iproc = proc_elastic;
        cout << " elastic "<<endl;
    }
    if (proc == "dis" || proc == "DIS") {
        iproc = proc_dis;
        cout << " deep inelastic "<<endl;
    }
    if (iproc == proc_undef) cout << " undefined ! (WARNING!)"<<endl;

    ebeam_central = eb;
    theta_central = theta;
    mass_tgt = masst;
    thmin = thmi;
    thmax = thma;
    phmin = phmi;
    phmax = phma;
    epmin = epmi;
    epmax = epma;

    did_init = kTRUE;

    return 1;
}

void hamcKine::SetDisDef(Float_t xlo, Float_t xhi, Float_t qslo, Float_t wslo) {
  xbjlo = xlo;
  xbjhi = xhi;
  qsqlo = qslo;
  wsqlo = wslo;

}

Int_t hamcKine::CheckInit() {
   if (!did_init) {
    cout <<"hamcKine::ERROR:  Did not initialize the class !"<<endl;
    return -1;
   }
   return 1;
}

Int_t hamcKine::Generate(hamcExpt *expt) {

// It's assumed the beam was already generated if it exists.

  beam = expt->event->beam;
  Float_t eb,E0;
  if (beam) {
    eb = beam->GetEnergy();  // beam energy this event
    E0 = beam->GetE0();
    // includes fluctuations and eloss before scattering.
  } else {
    eb = ebeam_central;
    E0 = ebeam_central;
  }

  hamcEloss *eloss = expt->physics->eloss;
  Float_t dE = 0;

// The internal is for t_equivalent, so it's
// what to subtract before and after scattering.

  if (eloss) dE = eloss->GetDeIntern() 
                + eloss->GetDeExternOut() 
                + eloss->GetDeIonOut();

// Add kick if we are iterating on energy
//   iteration = expt->iteration;
//   if (iteration == 1) {
//       dE = dE - dP0_iter*E0;
//   }
  
  if(Generate(eb, dE) == -1)  //no dis event found.
    return -1;
  if(energy <= (eprime+dE_after))      //physicslly not acceptable.
    return -1;

    return 1;

 

}

Int_t hamcKine::Generate(Float_t eb, Float_t dE) {

// To generate a track for an event.

  Clear();
  CheckInit();

  ebeam = eb;    energy = ebeam;
  dE_after = dE;

  if (iproc == proc_dis && energy*energy < wsqlo) return -1;

  if ( did_init == kFALSE ) {
    cout << "hamcKine::ERROR: uninitialized !"<<endl;
    return -1;
  }

// Generate scattering angle uniform in cos(theta)

  Float_t ctheta, cthmin, cthmax;

  cthmin = TMath::Cos(thmin);
  cthmax = TMath::Cos(thmax);

  ctheta = cthmin + (cthmax-cthmin)*gRandom->Rndm();

  theta = TMath::ACos(ctheta);

// Likewise, generate azimuthal angle (radians) in lab-frame
  phi = phmin + (phmax-phmin)*gRandom->Rndm();   

  if (iproc == proc_elastic) GenerateElastic();

  if (iproc == proc_dis) {
      if (GenerateDis() == -1) {
	//     cout << "hamcKine::WARNING: inf. loop? "<<ebeam<<endl;
        return -1;
      }
  }
  return 1;

}

Int_t hamcKine::GenerateElastic() {

  eprime = 0;

  if (mass_tgt == 0) return -1;

  // When iterating over incident angles, scattering angle is not simply theta, but the 
  //relative angle between incoming angle and outgoing angle;

  Float_t theta1, phi1, costheta;
  Float_t sts1, sps1, cts1,cps1, sts2, sps2, cts2, cps2;

  costheta = TMath::Cos(theta);
  theta1 = beam->dtheta_iter;
  phi1 = beam->dphi_iter;

  sts1 = TMath::Sin(theta1);
  sps1 = TMath::Sin(phi1);
  cts1 = TMath::Cos(theta1);
  cps1 = TMath::Cos(phi1);
  sts2 = TMath::Sin(theta);
  sps2 = TMath::Sin(phi);
  cts2 = TMath::Cos(theta);
  cps2 = TMath::Cos(phi);
  
  if ((theta1 ==0)&&(phi1 == 0))
    eprime = ebeam / ( 1 + ((ebeam/mass_tgt) * 
			    (1 - TMath::Cos(theta))) );
  else {
    
    if (!theta1)
      costheta = 1-((sts1 - sts2*cps2)*(sts1 - sts2*cps2)+(0-sts2*sps2)*(0-sts2*sps2)+(cts1-cts2)*(cts1-cts2))/2.;

    if (!phi1)
      costheta = 1-((0-sts2*cps2)*(0-sts2*cps2)+(sps1-sts2*sps2)*(sps1-sts2*sps2)+(cps1-cts2)*(cps1-cts2))/2.;

    eprime = ebeam / (1+((ebeam/mass_tgt)*(1-costheta)));
  }

  pprime = TMath::Sqrt(eprime*eprime - mass_electron*mass_electron);
  erecoil = ebeam + mass_tgt - eprime;
  
  ComputeKine();

// subtract radiative loss after scattering.
  eprime = eprime - dE_after;  

  return 1;
}

Int_t hamcKine::GenerateDis() {

// Loop until find an event within the cuts that define DIS

  Int_t maxloop = 100;
  Int_t iloop = 0;

  while (iloop++ < maxloop) {

    eprime = epmin + (epmax - epmin)*gRandom->Rndm();

    ComputeKine();

    eprime = eprime - dE_after;

 // Impose cuts that define DIS here
    if (x > xbjlo && x < xbjhi && qsq > qsqlo 
	&& wsq > wsqlo) return 1;     // Found a DIS event
  }

  return -1;

}

 
Int_t hamcKine::ComputeKine() {


  //  qsq = 2*ebeam*eprime*(1 - TMath::Cos(theta));

  Float_t px1,py1,pz1,px2,py2,pz2;
  Float_t sts, cts, sps, cps;

  px1 = beam->GetPx();
  py1 = beam->GetPy();
  pz1 = beam->GetPz();

  sts = TMath::Sin(theta);
  cts = TMath::Cos(theta);
  sps = TMath::Sin(phi);
  cps = TMath::Cos(phi);

  px2 = pprime*sps*sts;
  py2 = pprime*cps*sts;
  pz2 = pprime*cts;

  qsq = -1.0*((ebeam-eprime)*(ebeam-eprime)-((px1-px2)*(px1-px2)+(py1-py2)*(py1-py2)+(pz1-pz2)*(pz1-pz2)));

  Float_t mass = mass_tgt;
  if (iproc == proc_dis) mass = mass_proton;
  
  wsq = mass*mass + 2*mass*(ebeam-eprime) - qsq;
  x = qsq/2/mass/(ebeam-eprime);
  y = (ebeam - eprime) / ebeam;
  bigy = ( 1 - (1-y)*(1-y) ) / ( 1 + (1-y)*(1-y) );

  return 1;

}


void hamcKine::Print() {

  cout << "hamcKine: Print: "<<endl;
  cout << "energy "<<energy<<"   theta "<<theta<<"  phi "<<phi<<endl;
  cout << "qsq  "<<qsq<<"   wsq "<<wsq<<endl;

}


  
