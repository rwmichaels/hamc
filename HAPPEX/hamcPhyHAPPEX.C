//  hamcPhyHAPPEX   -- class for the physics of HAPPEX
//  R. Michaels  Nov 2008
//  Was extracted from gener_cone and modified a bit -- needs some checking.


#include "hamcPhyHAPPEX.h"
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
ClassImp(hamcPhyHAPPEX)
#endif


hamcPhyHAPPEX::hamcPhyHAPPEX() : hamcPhysics()
{
  phy_name = "HAPPEX physics";
  scatt_process = "elastic";
  do_radiate = kTRUE;
}


hamcPhyHAPPEX::~hamcPhyHAPPEX() { }


Int_t hamcPhyHAPPEX::Init(hamcExpt* expt) {

  didinit = kTRUE;

  hamcPhysics::Init(expt);


  return 1;
}


Int_t hamcPhyHAPPEX::Generate(hamcExpt *expt) {

   Float_t energy = expt->physics->kine->energy;
   Float_t theta = expt->physics->kine->theta;
   Float_t qsq = expt->physics->kine->qsq;

// Compute the cross section for HAPPEX
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



/* Returns the e-p elastic cross section in barn/sr. 
   E_beam in GeV, theta in deg, Q_sqr in GeV^2 

 (changed to barn/sr by R. Michaels, to be consistent with other generators) */

Float_t hamcPhyHAPPEX::sig_elas_H(Float_t E_beam, Float_t theta, Float_t Q_sqr)
{
  Float_t d_sig, tau, sin2, cos2, G_E2, G_M2;

  sin2 = pow(sin(theta/2.0),2);
  cos2 = pow(cos(theta/2.0),2);
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

 Float_t hamcPhyHAPPEX::GEp_dipol(Float_t Q_sqr)
{
  
     return 1. / pow( 1. + Q_sqr / MV2 , 2.);
}

/* The function GMp_dipol() returns the proton dipol 
   magnetic form factor as function of Q_sqr(in GeV).  */

Float_t hamcPhyHAPPEX::GMp_dipol(Float_t Q_sqr)
{
  
     return mu_p / pow( 1. + Q_sqr / MV2 , 2.);
}

/* The function GEn_dipol() returns the neutron dipol 
   electric form factor as function of Q_sqr(in GeV). 
   The numerical coefficient comes from the Galster
   parametrization */

Float_t hamcPhyHAPPEX::GEn_dipol(Float_t Q_sqr)
{
     Float_t tau;
  
     tau = Q_sqr/4/pow(Mp,2);
     return -1.* mu_n*tau / pow( 1.+ Q_sqr/MV2,2.) / (1+tau*lambda_n);
}

/* The function GMn_dipol() returns the neutron dipol 
   magnetic form factor as function of Q_sqr(in GeV).  */

Float_t hamcPhyHAPPEX::GMn_dipol(Float_t Q_sqr)
{
     return mu_n / pow( 1. + Q_sqr/MV2 , 2.);
}


Float_t hamcPhyHAPPEX::GAp3(Float_t Q_sqr)
{
     return GAp_dipol(Q_sqr)/2.;
}

Float_t hamcPhyHAPPEX::GAp8(Float_t Q_sqr)
{
    double M_A=1.001, F=0.463, D=0.804;

    return (3.*F-D)/2./sqrt(3.)/pow(1.+Q_sqr/pow(M_A, 2.),2.);
}

Float_t hamcPhyHAPPEX::GAp_dipol(Float_t Q_sqr)
{
     double M_A=1.001, gA=-1.2695;

     if(Q_sqr<0) Q_sqr = 0.0;

     return -gA/pow(1.+Q_sqr/pow(M_A, 2.),2.);
}



/***********************************************************************************/
Float_t hamcPhyHAPPEX::asym_H(Float_t theta, Float_t Q2)
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

Int_t hamcPhyHAPPEX::Drate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec) {

  tlen = tlen*100;   // need cm
  //cout<<"tlen="<<tlen<<", tdens = "<<tdens<<", anum="<<anum<<", crsec = "<<crsec<<endl;
  Float_t avg_omega = 0.0059;
  drate = 6.25e12 * crsec * 0.602 * tlen * tdens * avg_omega / anum; //Hz/uA
  return 1;
}
