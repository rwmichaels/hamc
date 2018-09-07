//  hamcPhyWater   -- class for the physics of Water Cell
//  R. Michaels  Sept 2018

#include "hamcPhyWater.h"
#include "hamcExpt.h"
#include "hamcTarget.h"
#include "hamcEvent.h"
#include "hamcBeam.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcInout.h"
#include "hamcKine.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h" 
#include "TROOT.h" 
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#ifndef NODICT
ClassImp(hamcPhyWater)
#endif


hamcPhyWater::hamcPhyWater() : hamcPhysics()
{
  phy_name = "Water physics";
  scatt_process = "elastic";
  do_radiate = kTRUE;
}


hamcPhyWater::~hamcPhyWater() { }


Int_t hamcPhyWater::Init(hamcExpt* expt) {

  didinit = kTRUE;

  hamcPhysics::Init(expt);

  Float_t theta_degr, theta_rad, ene, xcross, masstgt;
  Float_t frecoil, eprime, eprimeO16, eprimeH, deltaE, qsqloc;
  Float_t anum;

  if (quick_check) {

    // energies: PREX-I, PREX-II, CREX
    for (Int_t iiene=0; iiene<3; iiene++) {

      if(iiene==0) {
        ene=1.063;
      }
      if(iiene==1) {
        ene=0.95;
      }
      if(iiene==2) {
        ene=2.2;
      } 

      theta_degr = 5.0;
      theta_rad = 3.1415926*theta_degr/180.;

      // inuc: 0 = O16, 1 = H
      for (Int_t inuc=0; inuc<2; inuc++) {

	if (inuc == 0) {
	  masstgt = 15.0;
          anum = 16;
	} else {
	  masstgt = 0.938;
          anum = 1;
	}

        frecoil = 1 + (ene/masstgt)*(1-TMath::Cos(theta_rad));

        eprime = ene/frecoil;
        qsqloc = 2*ene*eprime*(1-TMath::Cos(theta_rad));
        Float_t qinvf = TMath::Sqrt(qsqloc)/0.197;


        if (inuc == 0) {
          eprimeO16 = eprime;
          xcross = O16CrossSection(ene, theta_degr);
	} else {
          eprimeH = eprime; 
          deltaE = eprimeH - eprimeO16;
          xcross = sig_elas_H(ene,theta_degr, qsqloc); 
	}
      string nucname;
      nucname="Oxygen";
      if(inuc==1) nucname="Hydrogen";
      cout << "\n\n******   "<<nucname<<"   "<<ene<<" GeV   "<<theta_degr << " degrees   q = "<<qinvf<<" fm^-1   ***************************"<<endl;
      if (inuc == 1) cout << "deltaE (O16 to H)  = "<<-1*deltaE*1000<<"  MeV "<<endl<<endl;;

      Float_t tdens = expt->target->GetMtlDensity(inuc);  // tgt density (g/cm^3)
      Float_t tlen = expt->target->GetMtlEffLen(inuc);  // tgt len (m)
      tlen = tlen*100;                        // need cm
      Float_t current = expt->event->beam->beam_current;  // microAmps (uA)
      current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec
      Float_t domega = 0.00293; // 1 HRS
      Float_t rate = current * xcross * 0.602 * tlen * tdens * domega / anum;

      cout << "A = "<<anum<<"   density = "<<tdens<<" g/cm^3   Eff. length = "<<tlen<< "cm    Current = "<<current<<" e-/sec"<<endl;

      cout << "cross sect. = "<<xcross<<" barns/str    solid ang = "<<domega<<" str    RATE = "<<rate<<" Hz/uA "<<endl;


      }
    }

    exit(0);

  }
}


Int_t hamcPhyWater::Generate(hamcExpt *expt) {

   Float_t energy = expt->physics->kine->energy;
   Float_t theta = expt->physics->kine->theta;
   Float_t qsq = expt->physics->kine->qsq;

// Compute the cross section for Water
   crsec =  sig_elas_H(energy, theta, qsq);

// Compute the PV asymmetry 
   asymmetry = asym_H(theta, qsq);

// Compute the differential rate
   Float_t anum = expt->target->GetAscatt();
   Int_t mtl_idx = expt->target->GetMtlIndex();
   Float_t tdens = expt->target->GetMtlDensity(mtl_idx);  // tgt density (g/cm^3)
   Float_t tlen = expt->target->GetMtlLen(mtl_idx);  // tgt len (m)

   Drate(anum, tdens, tlen, crsec); //Hz/uA

   return OK;
}

Float_t hamcPhyWater::O16CrossSection(Float_t energy, Float_t angle) {
// This finds the O16 cross section, and only for a few Q^2 values
//  If you want Hydrogen, use sig_elas_H

// energy in GeV
// angle in degrees
// result is in barn/seradians 

  Int_t ldebug=0;

  Float_t pi = 3.1415926; 

  Float_t qsqloc,qinvf;

  qsqloc = 2*energy*energy*(1-TMath::Cos(angle*pi/180.));

  qinvf = TMath::Sqrt(qsqloc)/0.197;

  Float_t calcrsec, mott, form_factor;

  /*Mott cross section for point-like scattering = 
    (alpha*Z* (hc/2pi)*cos(theta/2))^2 / 400E^2(sin(theta/2))^4 */
  
  Float_t halfangle_rad = (angle/2)*(pi/180);
  
  Float_t sin4 = pow(sin(halfangle_rad),4);
  Float_t znuc = 8;

  Float_t d_sig, tau, sin2, cos2, G_E2, G_M2;

  mott = pow(((znuc*0.197*cos(halfangle_rad))/137),2)/ (400*pow(energy,2)*sin4);

  // comparing to the "mott" from hydrogen calculation ... and they agreed.

  Float_t ang_rad = 3.14159*angle/180.;

  sin2 = pow(sin(ang_rad/2.0),2);
  cos2 = pow(cos(ang_rad/2.0),2);
  Float_t GEp = GEp_dipol(qsqloc);
  Float_t GMp = GMp_dipol(qsqloc);
  G_E2 = pow(GEp,2);
  G_M2 = pow(GMp,2); 
  tau = qsqloc/4/pow(Mp,2);

#ifdef THING1
  cout << "GEp, etc "<<sin2<<"  "<<cos2<<"  "<<GEp<<"   "<<GMp<<"   "<<G_E2<<"   "<<G_M2<<endl;
  cout << "d_sig part 1 "<<1e-6*pow(Alpha,2)/4./pow(energy,2)/pow(sin2,2)/(1.+2.*energy/Mp*sin2)<<endl;
  cout << "dbg1 "<<pow(Alpha,2)/4.<<endl;
  cout << "dbg2 "<<energy<<"  "<<pow(Alpha,2)/4./pow(energy,2)<<endl;
  cout << "dbg3 "<<pow(Alpha,2)/4./pow(energy,2)/pow(sin2,2)<<endl;
  cout << "dbg4 "<<1/(1.+2.*energy/Mp*sin2)<<endl;
#endif

  d_sig = pow(Alpha,2)/4./pow(energy,2)/pow(sin2,2)/(1.+2.*energy/Mp*sin2)
         *((G_E2+tau*G_M2)/(1.+tau)*cos2 + 2.*tau*G_M2*sin2);

  d_sig = d_sig * hbc2 * 1e-6;   // converts to barns/sr

  //  cout << "d_sig all "<<d_sig<<endl;
   //  cout << "Zsq * mott "<<znuc*znuc*d_sig<<endl;

  if(ldebug) cout << "Mott  "<<znuc<<"  "<<energy<<"   "<<angle<<"   "<<mott<<endl;

  // Extremely crude lookup of FF-squared.  We only need 3 values at the moment.

  form_factor = 0.645; // it is FF^2
  if (qinvf > 0.46 && qinvf < 0.9) form_factor = 0.575;
  if (qinvf >= 0.9) form_factor = 0.087;

  if (ldebug) cout << "q "<<qinvf<<"   ff^2 = "<<form_factor<<endl;

  calcrsec = mott*form_factor; //result is in barn/seradians 

  if (ldebug) cout << "CalcCross:  "<<angle<<"  "<<energy<<"   "<<qinvf<<"   "<<form_factor<<"   "<<
               mott<<"   "<<calcrsec<<endl;

  return calcrsec;

}


/* Returns the e-p elastic cross section in barn/sr. 
   E_beam in GeV, theta in deg, Q_sqr in GeV^2 

 (changed to barn/sr by R. Michaels, to be consistent with other generators) */

Float_t hamcPhyWater::sig_elas_H(Float_t E_beam, Float_t theta, Float_t Q_sqr)
{
  Float_t d_sig, tau, sin2, cos2, G_E2, G_M2;

  Float_t ang_rad = 3.14159*theta/180.;

  sin2 = pow(sin(ang_rad/2.0),2);
  cos2 = pow(cos(ang_rad/2.0),2);
  Float_t GEp = GEp_dipol(Q_sqr);
  Float_t GMp = GMp_dipol(Q_sqr);
  G_E2 = pow(GEp,2);
  G_M2 = pow(GMp,2); 
  tau = Q_sqr/4/pow(Mp,2);
  
  d_sig = pow(Alpha,2)/4./pow(E_beam,2)/pow(sin2,2)/(1.+2.*E_beam/Mp*sin2)
         *((G_E2+tau*G_M2)/(1.+tau)*cos2 + 2.*tau*G_M2*sin2);

  d_sig = d_sig * hbc2 * 1e-6;   // converts to barns/sr

  return d_sig;
}


/* The function GEp_dipol() returns the proton dipol 
   electric form factor as function of Q_sqr(in GeV).  */

 Float_t hamcPhyWater::GEp_dipol(Float_t Q_sqr)
{
  
     return 1. / pow( 1. + Q_sqr / MV2 , 2.);
}

/* The function GMp_dipol() returns the proton dipol 
   magnetic form factor as function of Q_sqr(in GeV).  */

Float_t hamcPhyWater::GMp_dipol(Float_t Q_sqr)
{
  
     return mu_p / pow( 1. + Q_sqr / MV2 , 2.);
}

/* The function GEn_dipol() returns the neutron dipol 
   electric form factor as function of Q_sqr(in GeV). 
   The numerical coefficient comes from the Galster
   parametrization */

Float_t hamcPhyWater::GEn_dipol(Float_t Q_sqr)
{
     Float_t tau;
  
     tau = Q_sqr/4/pow(Mp,2);
     return -1.* mu_n*tau / pow( 1.+ Q_sqr/MV2,2.) / (1+tau*lambda_n);
}

/* The function GMn_dipol() returns the neutron dipol 
   magnetic form factor as function of Q_sqr(in GeV).  */

Float_t hamcPhyWater::GMn_dipol(Float_t Q_sqr)
{
     return mu_n / pow( 1. + Q_sqr/MV2 , 2.);
}


Float_t hamcPhyWater::GAp3(Float_t Q_sqr)
{
     return GAp_dipol(Q_sqr)/2.;
}

Float_t hamcPhyWater::GAp8(Float_t Q_sqr)
{
    double M_A=1.001, F=0.463, D=0.804;

    return (3.*F-D)/2./sqrt(3.)/pow(1.+Q_sqr/pow(M_A, 2.),2.);
}

Float_t hamcPhyWater::GAp_dipol(Float_t Q_sqr)
{
     double M_A=1.001, gA=-1.2695;

     if(Q_sqr<0) Q_sqr = 0.0;

     return -gA/pow(1.+Q_sqr/pow(M_A, 2.),2.);
}



/***********************************************************************************/
Float_t hamcPhyWater::asym_H(Float_t theta, Float_t Q2)
/* Returns the parity violating asymetry in e-p elastic scattering (raw).
   rho and kappa factors account for electroweak corrections.
   No strange quark contribution.*/
/* R. Michaels.  The units are fraction, NOT ppm, to be compatible
   with other generators */
{
  Float_t tau, eps, epsp, GEp, GMp, GEn, GMn, GA3, GA8, asym;

  GEp = GEp_dipol(Q2);
  GMp = GMp_dipol(Q2);
  GEn = GEn_dipol(Q2);
  GMn = GMn_dipol(Q2);
  GA3 = GAp3(Q2);
  GA8 = GAp8(Q2);
  tau = Q2/4./pow(Mp,2);
  eps = 1./( 1.+2.*(1.+tau)*pow(tan(theta/2.),2) );
  epsp = sqrt(tau*(1.+tau)*(1.-eps*eps));

  asym = -1.*GFermi*Q2/(4.*sqrt(2.)*mypi*Alpha)*
          ( rhop*(1.-4.*kappap*sw2)-4.*lambda1u-2.*lambda1d
	  - (rhop+2.*lambda1u+4.*lambda1d)*(eps*GEp*GEn+tau*GMp*GMn)/(eps*pow(GEp,2)+tau*pow(GMp,2))
	  - epsp*GMp/(eps*pow(GEp,2)+tau*pow(GMp,2))
	        *( (rho*(1.-4.*kappa*sw2)-lambda2u+lambda2d)*(-2.)*GA3+(lambda2u+lambda2d)*2.*sqrt(3)*GA8 )
	  );

  return asym;
}

Int_t hamcPhyWater::Drate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec) {

  tlen = tlen*100;   // need cm
  //cout<<"tlen="<<tlen<<", tdens = "<<tdens<<", anum="<<anum<<", crsec = "<<crsec<<endl;
  Float_t avg_omega = 0.0059;
  drate = 6.25e12 * crsec * 0.602 * tlen * tdens * avg_omega / anum; //Hz/uA
  return 1;
}
