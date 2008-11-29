//  hamcTrack   -- A track
//  R. Michaels  June 2008

#include "hamcTrack.h"
#include "hamcTrans.h"
#include "hamcSpecHRS.h"
#include "hamcAperture.h"
#include "hamcExpt.h"
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

hamcTrack::hamcTrack(string pd) : origin(ITARGET),pid(pd),E0(0),P0(0),P0sigma(0),energy(0),theta(0),phi(0),did_init(kFALSE),inaccept(kFALSE)
{
  trktype="none";
  tvect = new hamcTransVect();
  tvect_orig = new hamcTransVect();
  InitMass();
}

hamcTrack::hamcTrack(string pd, Float_t ee, Float_t x, Float_t th, Float_t y, Float_t ph, Float_t dp) : origin(ITARGET),pid(pd),E0(ee),P0sigma(0),energy(ee),theta(0),phi(0),inaccept(kFALSE) {
  trktype="none";
  tvect = new hamcTransVect(x, th, y, ph, dp);   // gets modified by transport
  tvect_orig = new hamcTransVect(x, th, y, ph, dp);  // origin of track.
  InitMass();
}

hamcTrack::~hamcTrack() {
  delete tvect;
  delete tvect_orig;
}

Int_t hamcTrack::InitMass() {

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

  inaccept = aperture->CheckAccept(tvect->GetX(), tvect->GetY());

  return inaccept;

}

void hamcTrack::MultScatt(const hamcExpt *expt, Int_t where) {

  Float_t radlen;

  if (where == ITARGET) {  // only choice so far

    Float_t tlen =  expt->target->GetLength();  // "L" of target
    // zscat varies from -L/2 to +L/2
    Float_t zscat = expt->target->GetZScatt();
    Float_t x0rad = expt->target->GetRadLength();
    radlen = 0;
    if (x0rad != 0) radlen = ((tlen/2)-zscat)/x0rad;
    MultScatt(radlen, ITARGET);     

  }

}

void hamcTrack::MultScatt(const hamcAperture *aperture, Int_t where) {
// Mult scattering defined at each aperture.
// It might depend on X,Y location in acceptance (e.g.
// a Paulesque collimator with multiple materials).

  if (!aperture) return;  

  Float_t radlen = aperture->GetRadLen(tvect->GetX(), tvect->GetY());

  if (debug) cout << "Mult scatt "<<where<<"  "<<radlen<<endl;

  if (radlen <= 0) return;   // no mult scattering

  MultScatt(radlen, where);

}


void hamcTrack::MultScatt(Float_t radlen, Int_t where) {
// Applies multiple scattering to track parameters.
// Resets track origin to present location.
// The transport model must provide transport from this new origin.

  if (P0 == 0) {   
    cout << "hamcTrack::MultScatt: ERROR: P0 = 0 ?"<<endl;
    return;
  }

  Float_t theta_sigma = (0.0136/P0) * TMath::Sqrt(radlen);

  tvect->AddToTheta(theta_sigma * gRandom->Gaus());
  tvect->AddToPhi(theta_sigma * gRandom->Gaus());

  *tvect_orig = *tvect;
  origin = where;

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
