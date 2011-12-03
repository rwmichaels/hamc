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
#include "hamcSpecHRS.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>
#include "TH1F.h"
#include <math.h>

using namespace std;

#ifndef NODICT
ClassImp(hamcKine)
#endif


hamcKine::hamcKine(): did_init(kFALSE)
{

  use_eloss = 1;  // default is to use Eloss, but it can be turned off here.

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
  Float_t ebeam = expt->event->beam->GetEnergy();  // this may change later (Eloss)
// theta_central, converted to radians, only for single-arm.
  Float_t theta = PI * expt->GetSpectrom(0)->GetScattAngle() / 180;
  Float_t masst = expt->target->GetMass();

  iteration = 0;

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

   eprime_gen = new TH1F("eprime_gen", "eprime generated for PVDIS", 100, 0.9*emin, 1.1*emax);

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
  track = expt->event->trackout[0];
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
  dE_int = 0;

// The internal is for t_equivalent, so it's
// what to subtract before and after scattering.
// Note, if you use the "beam" energy it already
// has had eloss applied.

  if (eloss) dE_int = eloss->GetDeIntern(); 

  iteration = expt->iteration;  

  if(Generate(eb, dE_int) == -1) {  //no dis event found.
    //    cout<<"no dis event found"<<endl;
    return -1;
  }
  if(energy <= (eprime+dE_int)) {     //physicslly not acceptable.
    //   cout<<"physically unacceptable"<<endl;
    return -1;
  }
  return 1;

 
}

Int_t hamcKine::Generate(Float_t eb, Float_t dE) {

// To generate a track for an event.

  Clear();
  CheckInit();

  ebeam = eb;    
  energy = ebeam;

  // FIXME: the following is redundant.  
  // dE, which is passed as an argument of this method, is ALREADY dE_int, a member of the class
  dE_int = dE;

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
	//	cout << "hamcKine::WARNING: inf. loop? "<<ebeam<<endl;
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
  

  if (iteration == 0)

    eprime = ebeam / ( 1 + ((ebeam/mass_tgt) * 
		    	    (1 - TMath::Cos(theta))) );

  else if (iteration == 1) {

    if (theta1 > 1.E-8)
      costheta = 1-((sts1 - sts2*cps2)*(sts1 - sts2*cps2)+(0-sts2*sps2)*(0-sts2*sps2)+(cts1-cts2)*(cts1-cts2))/2.;

   else if (phi1 > 1.E-8)
      costheta = 1-((0-sts2*cps2)*(0-sts2*cps2)+(sps1-sts2*sps2)*(sps1-sts2*sps2)+(cps1-cts2)*(cps1-cts2))/2.;

   else costheta = TMath::Cos(theta);

    eprime = ebeam / (1+((ebeam/mass_tgt)*(1-costheta)));

  }
 
  else 
    cout<<"ERROR:hamcKine:   Generating Elastic error, iteration>1"<<endl;

  //  pprime = TMath::Sqrt(eprime*eprime - mass_electron*mass_electron);
  erecoil = ebeam + mass_tgt - eprime;
  
  ComputeKine();

// Subtract radiative loss after scattering.
// (the beam already had it's energy subtracted, and further subtraction will happen
//  on the weay out of the target in GenerateOut called by hamcEvent)

  if (use_eloss) eprime = eprime - dE_int;  

  return 1;
}

Int_t hamcKine::GenerateOut(hamcExpt *expt) {

  // Computes the qsq observed.

  Int_t ldebug = 0;

  Float_t px1,py1,pz1,px2,py2,pz2;
  Float_t dE;

  beam = expt->event->beam;
  track = expt->event->trackout[0];  

  hamcEloss *eloss = expt->physics->eloss;

  if (!beam || !track) {
    cout << "Error in hamcKine::GenerateOut "<<beam<<"  "<<track<<endl;
    return -1;
  }

// Subtract energy losses on the way out of the target

  dE = 0;
  if (eloss) dE = eloss->GetDeExternOut() 
                + eloss->GetDeIonOut();

  if (use_eloss) eprime = eprime - dE;

  track->UpdateFourMom(eprime);  // because of this, must use trackout[0] above

  ebeam = beam->GetEnergy();  
  px1 = beam->GetNoMsPx();  // Want to use the no-Mult-Scatt variable for true Qsq
  py1 = beam->GetNoMsPy();  // since scatt angle was generated relative to Z-beam.
  pz1 = beam->GetNoMsPz();

  px2 = track->GetPx();
  py2 = track->GetPy();
  pz2 = track->GetPz();


  Float_t bene = expt->event->beam->GetEnergy();  // microAmps (uA)
  Float_t theta = PI * expt->GetSpectrom(0)->GetScattAngle() / 180;
  Int_t which_hrs = expt->GetSpectrom(0)->which_spectrom;             // affects angle convention

// Sign for phi(Transport), i.e. horizontal angle:
// Use the sign convention for HRS, that it is +1/-1 for R/L HRS.
// However, if we have the Warm Septum, there is only a R-HRS transport, so we need to keep sign +1.

  Float_t xsign = 1.0;
  if (which_hrs == LEFTHRS) xsign = -1.0;  // sign convention, see above
  if (expt->GetSpectrom(0)->IsWarmSeptum()) xsign = 1.0;

  Float_t mcph1, mcth1, mcp1;
  Float_t mcph, mcth, mcp;

  // The following is the "observed" qsq.
  // The incoming beam is assumed to go along Z axis (that's what Podd assumes).
  // But there are Elosses in the beam.
  // The outgoing track has all mult scattering and all Elosses

  mcph1 = track->tvect->GetPhi();    // Transport phi (tangent of horizontal angle)
  mcth1 = track->tvect->GetTheta();  // Transport theta (tangent of vertical angle)
  mcp1  = track->GetPmom();          // track momentum

  // bene = beam energy
  // The following is the correct formula because it uses MS applied to Transport angles

  qsq_obs =  2*bene*mcp1*(1-((TMath::Cos(theta)+(xsign*(TMath::Sin(theta))*mcph1))/(TMath::Sqrt(1+mcth1*mcth1+mcph1*mcph1))));

  mcph = track->ph0;
  mcth = track->th0;
  mcp  = track->GetPmom();

  qsq_atrk = 2*bene*mcp*(1-((TMath::Cos(theta)+(xsign*(TMath::Sin(theta))*mcph))/(TMath::Sqrt(1+mcth*mcth+mcph*mcph))));

// In Podd, this would be qsq_atrk for left HRS (L) with xsign = -1.  (For R-HRS, xsign = +1).
//T->Draw("EK_L.Q2:2*(3.484)*(L.gold.p)*(1-((TMath::Cos(14.0*3.14159/180))+(xsign*(TMath::Sin(14.0*3.14159/180.))*L.gold.ph))/(TMath::Sqrt(1+L.gold.th*L.gold.th+L.gold.ph*L.gold.ph)))>>hqlc",ccut);

  qsqfr = -1;
  if (qsq != 0) {
    qsqfr = (qsq_obs-qsq)/qsq;
    //   cout << "Fractional "<<qsqfr <<endl;
  }

  if (ldebug) {
     cout << "\nObserved Qsq "<<endl;
     cout << "Beam "<<ebeam<<"  "<<px1<<"  "<<py1<<"  "<<pz1<<endl;
     cout << "Track "<<eprime<<"  "<<px2<<"  "<<py2<<"  "<<pz2<<endl;
     cout << "Qsq_obs  = "<<qsq_obs<<"   diff "<<qsq_obs-qsq<<endl;
     cout << "qsq_atrk "<<qsq_atrk<<"  qsq alone "<<qsq<<endl;
     cout << "Phi "<<mcph1 << "  "<<mcph<<endl;
     cout << "theta "<<mcth1<<"   "<<mcth<<endl;
     cout << "mcp " << mcp1<<"    "<<mcp<<endl;
     cout << "qsqfr = "<<qsqfr<<endl;
     if (TMath::Sqrt(qsqfr*qsqfr) > 2e-2) cout << "Large qsqfr !!"<<endl;
  }


  return 1;
}


Int_t hamcKine::GenerateDis() {

// Loop until find an event within the cuts that define DIS

  Int_t maxloop = 100;
  Int_t iloop = 0;

  while (iloop++ < maxloop) {

    eprime = epmin + (epmax - epmin)*gRandom->Rndm(1);

    ComputeKine();

 // Comment: probably don't want to subtract dE_int here, see comment above
    eprime = eprime - dE_int;

    //    eprime_gen->Fill(eprime);
    //    return 1;

 // Impose cuts that define DIS here
    if ((x > xbjlo) && (x < xbjhi) && (qsq > qsqlo) 
 	&& (wsq > wsqlo)) {
      eprime_gen->Fill(eprime);
       return 1;     // Found a DIS event
     }
//     else 
//       cout<<"trial failed x="<<x<<" qsq="<<qsq<<" wsq="<<wsq<<" eprime="<<eprime<<endl;
  }
  cout<<"beam energy:"<<energy<<"hamcKine::infinite loop while generating"<<endl;
  return -1;

}

 
Int_t hamcKine::ComputeKine() {

  Int_t ldebug = 0;

  Float_t px1,py1,pz1,px2,py2,pz2;
  Float_t sts, cts, sps, cps;
  Float_t sts1,cts1,sps1,cps1;
  Float_t theta1 = beam->dtheta_iter;
  Float_t phi1 = beam->dphi_iter;

  px1 = beam->GetNoMsPx();  // Want to use the no-Mult-Scatt variable for true Qsq
  py1 = beam->GetNoMsPy();  // since scatt angle was generated relative to Z-beam.
  pz1 = beam->GetNoMsPz();

  sts  = TMath::Sin(theta);
  cts  = TMath::Cos(theta);
  sps  = TMath::Sin(phi);
  cps  = TMath::Cos(phi);
  sts1 = TMath::Sin(theta1);
  cts1 = TMath::Cos(theta1);
  sps1 = TMath::Sin(phi1);
  cps1 = TMath::Cos(phi1);

  pprime = TMath::Sqrt(eprime*eprime - mass_electron*mass_electron);

  //  cout<<"pprime:"<<pprime<<endl;

// This is again the no-mult-scatt version of the outgoing momentum.

  px2 = pprime*cps*sts;
  py2 = pprime*sps*sts;
  pz2 = pprime*cts;

  // This is the "physical" qsq.
  // It has no mult scattering (see comment above) but the Eloss on
  // the beam was applied, yet no Eloss on the outgoing particle.

  qsq = -1.0*((ebeam-eprime)*(ebeam-eprime)-((px1-px2)*(px1-px2)+(py1-py2)*(py1-py2)+(pz1-pz2)*(pz1-pz2)));

  Float_t x1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  Float_t x2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  Float_t cloc = (px1*px2 + py1*py2 + pz1*pz2)/(x1*x2);

  if (ldebug) {
    cout << "\n\nTrue Qsq "<<endl;
    cout << "Beam "<<ebeam<<"  "<<px1<<"  "<<py1<<"  "<<pz1<<endl;
    cout << "eprime "<<eprime<<endl;
    cout << "Track "<<pprime<<"  "<<px2<<"  "<<py2<<"  "<<pz2<<endl;
    cout << "Qsq = "<<qsq<<endl;
    cout << "Cosines "<<cts<<"  "<<cloc<<"  "<<cloc-cts<<endl;
  }

  Float_t mass = mass_tgt;
  if (iproc == proc_dis) mass = mass_proton;
  
  wsq = mass*mass + 2*mass*(ebeam-eprime) - qsq;
  x = qsq/2/mass/(ebeam-eprime);
  y = (ebeam - eprime) / ebeam;
  bigy = ( 1 - (1-y)*(1-y) ) / ( 1 + (1-y)*(1-y) );


  if (iteration == 0)
    scat_ang = theta;

  else if (iteration == 1) {
    if (theta1 > 1.E-8) {
      Float_t   costheta = 1-((sts1 - sts*cps)*(sts1 - sts*cps)+(0-sts*sps)*(0-sts*sps)+(cts1-cts)*(cts1-cts))/2.;
      scat_ang = TMath::ACos(costheta);
      //   cout<<"theta="<<theta<<"scat_ang="<<scat_ang<<endl;
    }
    else if (phi1 > 1.E-8) {
      Float_t  costheta = 1-((0-sts*cps)*(0-sts*cps)+(sps1-sts*sps)*(sps1-sts*sps)+(cps1-cts)*(cps1-cts))/2.;
      scat_ang = TMath::ACos(costheta);
    }
    else scat_ang = theta;
  }

  return 1;

}


void hamcKine::Print() {

  cout << "hamcKine: Print: "<<endl;
  cout << "energy "<<energy<<"   theta "<<theta<<"  phi "<<phi<<endl;
  cout << "qsq  "<<qsq<<"   wsq "<<wsq<<endl;

}
