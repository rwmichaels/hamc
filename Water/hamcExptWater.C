//  hamcExptWater   -- Water cell measurements
//  R. Michaels  Dec 2018

#include "hamcExptWater.h"
#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcKine.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcBeam.h"
#include "hamcExpt.h"
#include "hamcSpecHRS.h"
#include "hamcTgtWater.h"
#include "hamcPhyWater.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptWater)
#endif

using namespace std;

hamcExptWater::hamcExptWater() : hamcSingles("Water")
{
  physics = new hamcPhyWater();
  target  = new hamcTgtWater();
}

hamcExptWater::~hamcExptWater() {
}

Int_t hamcExptWater::Init(string sfile) {

// Defaults in case the setup file is empty
  hamcSingles::SetP0(0.95);
  hamcSingles::SetTheta(5.0);

  hamcSingles::Init(sfile);

  // Defaults (can be over-ridden by input file)
  static_cast<hamcSpecHRS*>(spectrom[0])->UseCollimator();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseColdSeptum();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseMatrixTrans();

  spectrom[0]->Init(this);
  spectrom[0]->Print();


  return OK;
}

void hamcExptWater::EventAnalysis() {

// Fill some histograms for diagnostic purposes.

  Int_t lprint = 0;  // to turn on(1) or off(0) local print

  if (lprint) cout << "Into hamcExptWater:: Analysis "<<endl;

// Note, "event" and "physics" are public members of this class or it's parent

  // Qsq 
  Float_t qsq = physics->kine->qsq;

  // Angles at target
  Float_t th0 = event->trackout[0]->th0;
  Float_t ph0 = event->trackout[0]->ph0;

  // Transport variables at the focal plane
  Float_t x = event->trackout[0]->xtrans;
  Float_t th = event->trackout[0]->thtrans;
  Float_t y = event->trackout[0]->ytrans;
  Float_t ph = event->trackout[0]->phtrans;
  Float_t pmom = event->trackout[0]->GetPmom();

  Float_t theta_central = PI * GetSpectrom(0)->GetScattAngle() / 180;
  Int_t which_hrs = GetSpectrom(0)->which_spectrom;             // affects angle convention
  Float_t xsign = 1.0;
  if (which_hrs == LEFTHRS) xsign = -1.0;  // sign convention for L,R HRS.

  Float_t ebeam = event->beam->GetEnergy();  // energy of the beam right before scattering (contains Elosses)

  if (lprint) {
    cout << "\n\nAngles at target "<<th0<<"  "<<ph0<<endl;
    cout << "Transport vector at focal plane "<<x<<"  "<<y<<"  "<<th<<"  "<<ph<<endl;
    cout << "Qsq "<<qsq<<endl;

    cout << "inputs :  theta_central "<<theta_central<<"    Ebeam = "<<ebeam<<endl;
    cout << "angles at target  "<<ph0<<"  "<<th0<<"   Momentum "<<pmom<<endl;
    if (pmom > ebeam) cout << "Hmmm.  This should not happen ! "<<endl;
  }

  // You might want to make sure X,Y are extrapolated to the right plane.

  if (event->inaccept == 1) { // track is in focal plane acceptance

  }
  
}

