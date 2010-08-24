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

  qsq1 = new TH1F("qsq1","Qsq (weighted, in accept)",200,0.3,1.0);
  qsq2 = new TH1F("qsq2","Qsq (weighted, in accept, det accept cut)",200,0.3,1.0);
  qsq3 = new TH1F("qsq3","Qsq with cut1 ",200,0.3,1.0);
  qsq4 = new TH1F("qsq4","Qsq with cut1 ",200,0.3,1.0);
  hxy1 = new TH2F("hxy1","X vs Y in focal plane",100,-0.7,0.3,100,-0.06,0.06);
  hxy2 = new TH2F("hxy2","X vs Y in focal plane (det accept cut)",100,-0.7,0.3,100,-0.06,0.06);
  hxy3 = new TH2F("hxy3","X vs Y in focal plane, unweighted",100,-0.7,0.3,100,-0.06,0.06);
  hxy4 = new TH2F("hxy4","X vs Y in focal plane, det accept cut, unweighted",100,-0.7,0.3,100,-0.06,0.06);

  return OK;
}

void hamcExptHAPPEX::EventAnalysis() {

// Fill some histograms for diagnostic purposes.
// THIS IS JUST AN EXAMPLE for now.

  Int_t lprint = 0;  // to turn on(1) or off(0) local print

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

  // cuts to define detector location
  // (just an example; actually should make a trapezoid cut in the plane)

  // extracted below by plotting L.tr.x:L.tr.y
  // (x2,y2) = (-0.05,-0.6), (x1,y1)=(-0.01,0.25)
  Float_t xhi = 0.25, xlo = -0.6;
  Float_t y11 = -0.01, y12 = -0.05;
  Float_t y21 = 0.01, y22 = 0.05;

  // get the eqn of the line that bounds y
  // ylo: (-0.05,-0.01)
  Float_t slp1 = (y11-y12)/(xhi-xlo);
  Float_t ylo = y11 + slp1*(x-xhi);

  // yhi: (0.01, 0.05)
  Float_t slp2 = (y21-y22)/(xhi-xlo);
  Float_t yhi = y21 + slp2*(x-xhi);

//   Float_t ylo = -0.05;    
//   Float_t yhi =  0.05;    
  
  // Cross section for weight factor
  Float_t crsec = physics->GetCrossSection();
  Float_t wt = 1e5*crsec;

  if (lprint) {
    cout << "Angles at target "<<th0<<"  "<<ph0<<endl;
    cout << "Transport vector at focal plane "<<x<<"  "<<y<<"  "<<th<<"  "<<ph<<endl;
    cout << "Qsq "<<qsq<<endl;
    cout << "Cross section "<<crsec<<endl<<endl;
  }

  // You might want to make sure X,Y are extrapolated to the right plane.

  if (event->inaccept == 1) { // track is in focal plane acceptance

    hxy1->Fill(x,y,crsec);
    qsq1->Fill(qsq,wt);
    
    hxy3->Fill(x,y);

    if (x > xlo && x < xhi && y > ylo && y < yhi) {  // track in detector
      
      hxy2->Fill(x,y,crsec);
      qsq2->Fill(qsq,wt);

      hxy4->Fill(x,y);
    }
  }
  
}

