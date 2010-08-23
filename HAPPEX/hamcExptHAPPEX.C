//  hamcExptHAPPEX   -- HAPPEX Experiment.
//  R. Michaels  Dec 2008

#include "hamcExptHAPPEX.h"
#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcKine.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcBeam.h"
#include "hamcExpt.h"
#include "hamcSpecHRS.h"
#include "hamcTgtHAPPEX.h"
#include "hamcPhyHAPPEX.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptHAPPEX)
#endif

using namespace std;

hamcExptHAPPEX::hamcExptHAPPEX() : hamcSingles("HAPPEX")
{
  physics = new hamcPhyHAPPEX();
  target  = new hamcTgtHAPPEX();
}

hamcExptHAPPEX::~hamcExptHAPPEX() {
}

Int_t hamcExptHAPPEX::Init(string sfile) {

// Defaults in case the setup file is empty
  hamcSingles::SetP0(3.176);
  hamcSingles::SetTheta(6.106);  // for R-HRS

  hamcSingles::Init(sfile);

  // Defaults (can be over-ridden by input file)
  static_cast<hamcSpecHRS*>(spectrom[0])->UseCollimator();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseColdSeptum();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseMatrixTrans();

  spectrom[0]->Init(this);
  spectrom[0]->Print();

  qsq1 = new TH1F("qsq1","Qsq with cut1 ",200,0.3,1.0);
  qsq2 = new TH1F("qsq2","Qsq with cut1 ",200,0.3,1.0);
  qsq3 = new TH1F("qsq3","Qsq with cut1 ",200,0.3,1.0);
  qsq4 = new TH1F("qsq4","Qsq with cut1 ",200,0.3,1.0);
  hxy1 = new TH2F("hxy1","X vs Y in focal plane",100,-1,0.5,100,-0.2,0.2);
  hxy2 = new TH2F("hxy2","X vs Y in focal plane (cut)",100,-1,0.5,100,-0.2,0.2);

  return OK;
}

void hamcExptHAPPEX::EventAnalysis() {

// Fill some histograms for diagnostic purposes.
// THIS IS JUST AN EXAMPLE for now.

  Int_t lprint = 0;  // to turn on(1) or off(0) local print

// cuts to define detector location
// (just an example; actually should make a trapezoid cut in the plane)
  Float_t xlo = -0.6;    
  Float_t xhi =  0.25;    
  Float_t ylo = -0.05;    
  Float_t yhi =  0.05;    
 
  if (lprint) cout << "Into hamcExptHAPPEX:: Analysis "<<endl;

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

  // Cross section for weight factor
  Float_t crsec = physics->GetCrossSection();

  if (lprint) {
    cout << "Angles at target "<<th0<<"  "<<ph0<<endl;
    cout << "Transport vector at focal plane "<<x<<"  "<<y<<"  "<<th<<"  "<<ph<<endl;
    cout << "Qsq "<<qsq<<endl;
    cout << "Cross section "<<crsec<<endl<<endl;
  }

  // You might want to make sure X,Y are extrapolated to the right plane.

  if (event->inaccept == 1) { // track is in focal plane acceptance

     hxy1->Fill(x,y,crsec);
     qsq1->Fill(qsq,crsec);

     if (x > xlo && x < xhi && y > ylo && y < yhi) {  // track in detector

       hxy2->Fill(x,y,crsec);
       qsq2->Fill(qsq,crsec);

     }
  }


}

