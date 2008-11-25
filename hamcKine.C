//  hamcKine   -- Kinematics generator.
//  Member of hamcPhysics
//  Kinematics are generated uniform in phase space.
//  Works for elastic and DIS.

//  R. Michaels  Nov 2008


#include "hamcKine.h"
#include "hamcExpt.h"
#include "hamcPhysics.h"
#include "hamcBeam.h"
#include "hamcRad.h"
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

  dP0_iter = 0;

  parser.Load(expt->inout->GetStrVect("kick:track"));
  parser.Print();
  if (parser.IsFound("P0")) {
    if (expt->inout->numiter > 1) dP0_iter = parser.GetData();  
    cout << "Will iterate P0 by "<<100*dP0_iter<<" %"<<endl;
  }   


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

    SetThetaTable();

    did_init = kTRUE;

    return 1;
}

void hamcKine::SetDisDef(Float_t xlo, Float_t xhi, Float_t qslo, Float_t wslo) {
  xbjlo = xlo;
  xbjhi = xhi;
  qsqlo = qslo;
  wsqlo = wslo;

}

void hamcKine::SetThetaTable() {
// Setup table of theta values such that 
// solid angle sin(theta)*dtheta*dphi
// is uniformly filled

  Float_t th,x;
  Int_t i,k,num;

  Float_t cdiv = 200;
  Int_t MAX = (Int_t)(MAXCELL/cdiv);

  numtcell = 0;

  for (i=0; i<MAX; i++) {

    x = (Float_t)i/((Float_t)MAX);
    th = thmin + x * (thmax-thmin);
    thetacell.push_back(th);

    num = ((Int_t)(cdiv * TMath::Sin(th)));

    for (k=0; k<num; k++) tcellnum.push_back(i);

    numtcell += num;
 
    if (numtcell > MAXCELL) {
      // This should never happen.  If it does, ask me.
      cout << "hamcKine:: ERROR:  trying to create too many cells."<<endl;
      exit(0);
    }

  }

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

  hamcBeam *beam = expt->event->beam;
  Float_t eb,E0;
  if (beam) {
    eb = beam->GetEnergy();  // beam energy this event
    E0 = beam->GetE0();
    // includes fluctuations and eloss before scattering.
  } else {
    eb = ebeam_central;
    E0 = ebeam_central;
  }
  hamcRad *rad = expt->physics->radiation;
  Float_t dE = 0;
  if (rad) dE = rad->GetDeIntern() + rad->GetDeExternOut(); 

// Add kick if we are iterating on energy
  if (expt->iteration == 1) {
      dE = dE - dP0_iter*E0;
  }
  
  return Generate(eb, dE);

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

// Generate scattering angle (radians) in lab-frame.
 
  theta = 0.0001;
  if (tcellnum.size()>0 && thetacell.size()>0) {
    Float_t x = ((Float_t)numtcell)*gRandom->Rndm();
    Int_t icell = (Int_t) x;
    if (icell >= 0 && icell < thetacell.size()) {
        Int_t index = tcellnum[icell];
        if (index >= 0 && index < thetacell.size()) {  
            theta = thetacell[index];
	}
    }
  } else {
    cout << "hamcTrackOut::ERROR in theta generation"<<endl;
  }

// Likewise, generate azimuthal angle (radians) in lab-frame
  phi = phmin + (phmax-phmin)*gRandom->Rndm();   

  if (iproc == proc_elastic) GenerateElastic();

  if (iproc == proc_dis) {
      if (GenerateDis() == -1) {
	//        cout << "hamcKine::WARNING: inf. loop? "<<ebeam<<endl;
        return -1;
      }
  }
  return 1;

}

Int_t hamcKine::GenerateElastic() {

  eprime = 0;

  if (mass_tgt == 0) return -1;

  eprime = ebeam / ( 1 + ((ebeam/mass_tgt) * 
		         (1 - TMath::Cos(theta))) );

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

 // Impose cuts that define DIS here
    if (x > xbjlo && x < xbjhi && qsq > qsqlo 
	&& wsq > wsqlo) return 1;     // Found a DIS event

  }

  return -1;

}

 
Int_t hamcKine::ComputeKine() {

  qsq = 2*ebeam*eprime*(1 - TMath::Cos(theta));

  Float_t mass = mass_tgt;
  if (iproc == proc_dis) mass = mass_proton;
  
  wsq = 2*mass*(ebeam-eprime) - qsq;
  x = qsq/2/mass/(ebeam-eprime);
  y = (ebeam - eprime) / ebeam;
  bigy = ( 1 - (1-y)*(1-y) ) / ( 1 + (1-y)*(1-y) );

  return 1;

}
  

  
