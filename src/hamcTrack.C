//  hamcTrack   -- A track
//  R. Michaels  Feb 2009

#include "hamcMultScatt.h"
#include "hamcTrack.h"
#include "hamcTrans.h"
#include "hamcSpecHRS.h"
#include "hamcAperture.h"
#include "hamcExpt.h"
#include "hamcPhysics.h"
#include "hamcEloss.h"
#include "hamcTarget.h"
#include "hamcBeam.h"
#include "TRandom.h"
#include "TMath.h"
#include "Rtypes.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTrack)
#endif

hamcTrack::hamcTrack(string pd) : origin(ITARGET),pid(pd),E0(0),P0(0),P0sigma(0),energy(0),theta(0),phi(0),ms_collim(0),did_init(kFALSE),inaccept(kFALSE)
{
  trktype="none";
  tvect = new hamcTransVect();
  tvect_orig = new hamcTransVect();
  InitVars();
}

hamcTrack::hamcTrack(string pd, Float_t ee, Float_t x, Float_t th, Float_t y, Float_t ph, Float_t dp) : origin(ITARGET),pid(pd),E0(ee),P0sigma(0),energy(ee),theta(0),phi(0),ms_collim(0),inaccept(kFALSE) {
  trktype="none";
  tvect = new hamcTransVect(x, th, y, ph, dp);   // gets modified by transport
  tvect_orig = new hamcTransVect(x, th, y, ph, dp);  // origin of track.
  InitVars();
}

hamcTrack::~hamcTrack() {
  delete tvect;
  delete tvect_orig;
}

hamcTrack& hamcTrack::operator=(const hamcTrack& rhs) {
  *tvect = *rhs.tvect;
  *tvect_orig = *rhs.tvect_orig;
  xtrans = rhs.xtrans;
  thtrans = rhs.thtrans;
  ytrans = rhs.ytrans;
  phtrans = rhs.phtrans;
  dpptrans = rhs.dpptrans;
  ztrans = rhs.ztrans;
  x0 = rhs.x0;
  th0 = rhs.th0;
  y0 = rhs.y0;
  ph0 = rhs.ph0;
  dpp0 = rhs.dpp0;
  z0 = rhs.z0;
  th_ms = rhs.th_ms;
  ph_ms = rhs.ph_ms;
  thtgt = rhs.thtgt;
  phtgt = rhs.phtgt;
  E0 = rhs.E0;
  P0 = rhs.P0;
  E0sigma = rhs.E0sigma;
  P0sigma = rhs.P0sigma;
  energy = rhs.energy;
  pmom = rhs.pmom;
  mass = rhs.mass;
  theta = rhs.theta;
  phi = rhs.phi;
  theta_deg = rhs.theta_deg;
  phi_deg = rhs.phi_deg;
  dP0_iter = rhs.dP0_iter;
  pnoms_x = rhs.pnoms_x;
  pnoms_y = rhs.pnoms_y;
  pnoms_z = rhs.pnoms_z;
  plab_x = rhs.plab_x;
  plab_y = rhs.plab_y;
  plab_z = rhs.plab_z;
  xdet = rhs.xdet;
  ydet = rhs.ydet;
  return *this;
}


Int_t hamcTrack::InitVars() {

  use_mscat = 1;  // default should be 1, but can turn it off here.

  mass=0;
  if (pid=="electron") mass=5.11e-4;
  if (pid=="pion") mass=0.1396;
  if (pid=="proton") mass=0.938;

  if (mass == 0) {
    cout << "WARNING: hamcTrack: mass = 0  ???"<<endl;
    return ERROR;
  }

  return OK;
}

Int_t hamcTrack::Eloss(const hamcExpt *expt, const hamcAperture *aperture, Int_t where) {

 // For an energy loss in a location other than target.
 // The target radiation is handled in event generation.

  Float_t radlen = aperture->GetRadLen(tvect->GetX(), tvect->GetY());

  if (radlen <= 0) return 0;   // not in the A_T hole

  Float_t dE = expt->physics->eloss->GetDeDx(radlen, where);

  if (P0 != 0) {
    Float_t P = P0*(1+tvect->GetDpp());
    P = P - dE;  // assumes rel. electron
    Float_t dpn = (P-P0)/P0;
    tvect->PutDpp(dpn);
    *tvect_orig = *tvect;   // reset the original tvect.
  }

  return OK;

}

Int_t hamcTrack::UpdateFourMom(Float_t newE) {

// Update the four-momentum given an overall energy change.

  Float_t x1,x2,escale;

  x1 = TMath::Sqrt(energy*energy - mass*mass);

// These are cases (~0.2 %) where the Eloss is large but also slightly
// wrong since Eloss assumed initial energy = E0 but it might have been
// less than E0 if it was after another Eloss or if after scattering.
// Doesn't matter much; these are rejected by momentum acceptance.
// To simply reject these events we set newE, pmom, and dP/P here.

  if (newE < mass || x1 == 0) {
    newE = 0.001;
    pmom = 0.001;
    tvect->PutDpp(-0.999);   // lost nearly all energy --> reject
    *tvect_orig = *tvect;   // reset the original tvect.
    return -1; 
  }

  x2 = TMath::Sqrt(newE*newE - mass*mass);

  escale = x2/x1;

  energy = newE;
  plab_x = escale * plab_x;
  plab_y = escale * plab_y;
  plab_z = escale * plab_z;

  pmom = TMath::Sqrt(plab_x*plab_x + plab_y*plab_y + plab_z*plab_z);

  if (P0 != 0) {
    Float_t dpn = (pmom-P0)/P0;
    tvect->PutDpp(dpn);
    *tvect_orig = *tvect;   // reset the original tvect.
  }

  return 1;

}


Int_t hamcTrack::Init() {

  E0 = TMath::Sqrt(P0*P0 + mass*mass);
  energy = E0;  

  return OK;

} 

void hamcTrack::Print() {

  cout << "track type = "<<trktype<< "   pid = "<<pid<<endl;
  cout << "mass = "<<mass<<"   central energy = "<<E0<<"  GeV"<<endl;
  cout << "energy = "<<energy<<"    momentum = "<<pmom<<endl;
  tvect->Print();
  cout << "polar angles "<<theta<<"  "<<phi<<endl;  

}   

Int_t hamcTrack::Transport(const hamcTrans *trans, Int_t where) {

// Transport uses the hamcTrans object to transform this track
// from it's origin tvect_orig to "where".

     if (!trans) return ERROR;

     if (debug) cout << "\n\n -------------- \n Transform to "<<where<<endl;

     trans->TransForm(this, where);

     UpdateTrans();

     if (debug) {
         cout << "\n orig tvect "<<endl;
         tvect_orig->Print();
         cout << "\n tvect "<<endl;
         tvect->Print();
     }

     return OK;

}

Bool_t hamcTrack::InAccept(const hamcAperture *aperture) {
	
   if (!aperture) {   // no acceptance cut
    inaccept = kTRUE;
    return OK;  
  }

  inaccept = aperture->CheckAccept(this);  
 
  return inaccept;

}

void hamcTrack::MultScatt(const hamcExpt *expt, Int_t where) {

  Float_t radlen;

  if (where == ITARGET) {  

      /*
       // SPR 8/29:  I think this is antiquated
       // zscat varies from -L/2 to +L/2 + offset now
       
    Float_t tlen =  expt->target->GetLength();  // "L" of target
    // zscat varies from -L/2 to +L/2
    Float_t zscat = expt->target->GetZScatt();
    Float_t x0rad = expt->target->GetRadLength();
    radlen = 0;
    if (x0rad != 0) radlen = (((tlen/2)-zscat)/tlen) * x0rad;

    MultScatt(radlen, where);     
    */

      MultScatt(expt->target->GetPartialScatt(P0), where );

  }

  if (where == ITARGET_FULL) {  

      /*
    Float_t x0rad = expt->target->GetRadLength();
    MultScatt(x0rad, where);     
    */

      MultScatt(expt->target->GetFullScatt(P0), where );

  }

}

void hamcTrack::MultScatt(const hamcAperture *aperture, Int_t where) {
// Mult scattering defined at each aperture.
// It might depend on X,Y location in acceptance (e.g.
// a Paulesque collimator with multiple materials).

  if (!aperture) return;  


  //Float_t radlen = aperture->GetRadLen(tvect->GetX(), tvect->GetY());

  Double_t a = aperture->GetA(tvect->GetX(), tvect->GetY());
  Double_t z = aperture->GetZ(tvect->GetX(), tvect->GetY());
  Double_t t = aperture->Gett(tvect->GetX(), tvect->GetY());


  if (t <= 0){
      return;   // no mult scattering
  }

  hamcMultScatt *appscat = new hamcMultScatt(P0, t, a, z);

//  MultScatt(radlen, where);
  MultScatt(appscat, where);

  delete appscat;
}


void hamcTrack::MultScatt(Float_t radlen, Int_t where) {
    fprintf(stderr, "%s %s line %d:  ERROR:  Shouldn't be using this multiple scattering code anymore...\n",  __FILE__, __FUNCTION__, __LINE__ );
    exit(1);
// Applies multiple scattering to track parameters.
// Resets track origin to present location.
// The transport model must provide transport from this new origin.

  Int_t use_resol = 0;

  Float_t vresol,hresol,vkick,hkick;

  vresol = 0.002;        // resolution of the HRS in vertical angle
  hresol = 0.0006;       // ditton, horizontal

  if (P0 == 0) {   
    cout << "hamcTrack::MultScatt: ERROR: P0 = 0 ?"<<endl;
    return;
  }

  Float_t theta_sigma = (0.0136/P0) * TMath::Sqrt(radlen);

  Float_t dtheta1, dtheta2, prob, probsign, xsign;

  xsign = 1;

  prob = gRandom->Rndm();
  if (prob < 0.02) {
    probsign = gRandom->Rndm();
    if (probsign < 0.5) xsign = -1;
    dtheta1 = xsign * 5 * theta_sigma * gRandom->Rndm();  // flat tail
  } else {
    dtheta1 = theta_sigma * gRandom->Gaus();
  }

  if (use_mscat) tvect->AddToTheta(dtheta1);

  prob = gRandom->Rndm();
  if (prob < 0.02) {
    probsign = gRandom->Rndm();
    if (probsign < 0.5) xsign = -1;
    dtheta2 = xsign * 5 * theta_sigma * gRandom->Rndm();  // flat tail
  } else {
    dtheta2 = theta_sigma * gRandom->Gaus();
  }

  if (use_mscat) tvect->AddToPhi(dtheta2);

  *tvect_orig = *tvect;

  if (where == ITARGET_FULL) origin = ITARGET;

  if (where == ICOLLIM2) ms_collim=1;

  vkick = 0; 
  hkick = 0;
  if (where == ITARGET || where == ITARGET_FULL) {
     thtgt = tvect->GetTheta();
     phtgt = tvect->GetPhi();
    // Resolution smearing -- affects thtgt and phtgt variables.
    if (use_resol) {
       vkick = vresol*gRandom->Gaus();
       hkick = hresol*gRandom->Gaus();
       thtgt = thtgt + vkick;
       phtgt = phtgt + hkick;
     }
  }

  th_ms = tvect->GetTheta();
  ph_ms = tvect->GetPhi();

  plab_x = pnoms_x + pmom * (dtheta1 + hkick);
  plab_y = pnoms_y + pmom * (dtheta2 + vkick) ;
  plab_z = TMath::Sqrt(pmom*pmom - plab_x*plab_x - plab_y*plab_y);

}


void hamcTrack::UpdateTrans() {

  xtrans   = tvect->GetX();
  thtrans  = tvect->GetTheta();
  ytrans   = tvect->GetY();
  phtrans  = tvect->GetPhi();
  dpptrans = tvect->GetDpp();
  ztrans   = tvect->GetZ();

  xdet = xtrans/0.707;
  ydet = ytrans;


}

void hamcTrack::MultScatt(hamcMultScatt *ms, Int_t where) {
// Applies multiple scattering to track parameters.
// Resets track origin to present location.
// The transport model must provide transport from this new origin.

  Int_t use_resol = 0;

  Float_t vresol,hresol,vkick,hkick;

  vresol = 0.002;        // resolution of the HRS in vertical angle
  hresol = 0.0006;       // ditton, horizontal


  Float_t dtheta1, dtheta2;

  dtheta1 = ms->GenerateMSPlane();
  dtheta2 = ms->GenerateMSPlane();

  //  cout << "dtheta "<<dtheta1<<"  "<<dtheta2<<endl;

  if (use_mscat) tvect->AddToTheta(dtheta1);

  if (use_mscat) tvect->AddToPhi(dtheta2);

  *tvect_orig = *tvect;

  if (where == ITARGET_FULL) origin = ITARGET;

  if (where == ICOLLIM2) ms_collim=1;


  vkick = 0; 
  hkick = 0;
  if (where == ITARGET || where == ITARGET_FULL) {
     thtgt = tvect->GetTheta();
     phtgt = tvect->GetPhi();
    // Resolution smearing -- affects thtgt and phtgt variables.
    if (use_resol) {
       vkick = vresol*gRandom->Gaus();
       hkick = hresol*gRandom->Gaus();
       thtgt = thtgt + vkick;
       phtgt = phtgt + hkick;
     }
  }

  th_ms = tvect->GetTheta();
  ph_ms = tvect->GetPhi();

  // Following is obsolete.  MS should NOT be applied this way 
  // (plab is in different frame from Transport vector !)
  // Fortunately it does not affect qsq_obs anymore.

  plab_x = pnoms_x + pmom * (dtheta1 + hkick);
  plab_y = pnoms_y + pmom * (dtheta2 + vkick) ;
  plab_z = TMath::Sqrt(pmom*pmom - plab_x*plab_x - plab_y*plab_y);


}
