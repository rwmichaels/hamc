//  hamcExptPVDIS   -- PVDIS
//  D. Wang   R. Michaels  Jan. 2009

#include "hamcExptPVDIS.h"
#include "hamcExpt.h"
#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcKine.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcBeam.h"
#include "hamcSpecHRS.h"
#include "hamcAccAvg.h"
#include "hamcTgtPVDIS.h"
#include "hamcPhyPVDIS.h"
#include "hamcEloss.h"
#include "hamcInout.h"
#include "THaString.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptPVDIS)
#endif

using namespace std;

hamcExptPVDIS::hamcExptPVDIS() : hamcSingles("PVDIS")
{
  physics = new hamcPhyPVDIS();
  target  = new hamcTgtPVDIS();
}

hamcExptPVDIS::~hamcExptPVDIS() {
}

Int_t hamcExptPVDIS::Init(string sfile) {

// Defaults in case the setup file is empty
  hamcSingles::SetP0(1.04);
  hamcSingles::SetTheta(5.0);

  hamcSingles::Init(sfile);

  hpvd1 = new TH1F("hpvd1","Energies PVDIS",100,4,6.1);
  hpvd2 = new TH1F("hpvd2","Qsq PVDIS",100,0.5,2.2);
  hpvd3 = new TH1F("hpvd3","x PVDIS",100,0,1);
  hpvd4 = new TH2F("hpvd4","x,y PVDIS",100,-1,1,100,-0.3,0.3);

  hpvd5 = new TH1F("hpvd5","(accepted) Energies PVDIS",100,4,6.1);
  hpvd6 = new TH1F("hpvd6","(accepted) Qsq PVDIS",100,0.5,2.2);
  hpvd7 = new TH1F("hpvd7","(accepted) x PVDIS",100,0,1);
  hpvd8 = new TH2F("hpvd8","(accepted) x,y PVDIS",100,-1,1,100,-0.3,0.3);

  hpvd9 = new TH1F("hpvd9","Energy loss (in acceptance)",400,-0.02,2);
  hpvd10 = new TH1F("hpvd10","Depolar noscreen, 1",2000,-0.002,0.075);
  hpvd11 = new TH1F("hpvd11","Depolar noscreen 2",2000,-0.002,1);
  hpvd12 = new TH1F("hpvd12","(absolute) fractional energy loss ",100,-0.002,0.1);
  hpvd13 = new TH1F("hpvd13","Depolar screen, 1",2000,-0.002,0.075);
  hpvd14 = new TH1F("hpvd14","Depolar screen 2",2000,-0.002,1);
  hpvd15 = new TH1F("hpvd15","Diff: screen vs no screen",1000,-0.007,0.007);


  // Defaults (can be over-ridden by input file)
  static_cast<hamcSpecHRS*>(spectrom[0])->UseCollimator();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseWarmSeptum();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseMatrixTrans();

  spectrom[0]->Init(this);
  spectrom[0]->Print();

  return OK;
}

void hamcExptPVDIS::EventAnalysis() {

  hamcSingles::EventAnalysis();

  Int_t lprint=0;

  Float_t qsq = physics->kine->qsq;

  Float_t energy = physics->kine->energy;
  Float_t eprime = physics->kine->eprime;
  Float_t xval = physics->kine->x;
  Float_t p0 = GetSpectrom(0)->GetP0();

  Float_t theta = (180./3.14159)*physics->kine->theta;
  Float_t sth = TMath::Sin(physics->kine->theta);
  Float_t thwt = 1;
  if (sth > 0) thwt = 1./sth;

  Float_t th0 = event->trackout[0]->th0;
  Float_t ph0 = event->trackout[0]->ph0;

  Float_t x = event->trackout[0]->xtrans;
  Float_t th = event->trackout[0]->thtrans;
  Float_t y = event->trackout[0]->ytrans;
  Float_t ph = event->trackout[0]->phtrans;
  Float_t wt = 1e5*physics->GetCrossSection();
  Float_t th_deg = (180.0/PI) * event->trackout[0]->GetTheta();
  Float_t dExin = physics->eloss->GetDeExternIn();
  Float_t dExout = physics->eloss->GetDeExternOut();
  Float_t dEint = physics->eloss->GetDeIntern();
  Float_t dEtot = dExin + dExout + dEint;
  Float_t dE = dExin + 0.5*dEint;

  Float_t beampol = 0.89;

  // This is the formula 9.3 of Olsen & Maximom; no screening.
  Float_t depolnos = (dE*dE)*(1-(0.333*beampol*beampol))/(energy*energy + eprime*eprime - 0.666*energy*eprime);

  // Formulat 9.11 with screening
  Float_t psi1=20.8;
  Float_t psi2=20.2;

  Float_t depol = (dE*dE)*(psi1-(beampol*beampol*(psi1-0.666*psi2)))/(((energy*energy + eprime*eprime)*psi1) - 0.666*energy*eprime*psi2);

  Float_t fracE = (eprime - p0)/p0;

  hpvd1->Fill(energy);
  hpvd2->Fill(qsq);
  hpvd3->Fill(xval);
  hpvd4->Fill(x,y);

  if (fracE < 0) fracE = -1*fracE;
  Float_t wei = 6e-5*inout->GetWeight();

  if (event->inaccept == 1 && wei > 0) {

    if (lprint) {
      cout << "PVDIS event analysis " << energy <<"  "<<eprime<<"  "<<qsq<<"  "<<x<<"  "<<th_deg<<endl;
      cout << "Elosses "<<dExin <<"  "<<dExout<<"  "<<dEint<<"  "<<dEtot<<endl;
      cout << "Depol no screen "<<depolnos<<"    screen "<<depol<<endl;
      cout << "P0 "<<p0<<"   "<<fracE<<endl;
      cout << "Weight "<<wei<<endl;
    }

    hpvd5->Fill(energy, wei);
    hpvd6->Fill(qsq, wei);
    hpvd7->Fill(xval, wei);
    hpvd8->Fill(x,y);
    hpvd9->Fill(dEtot, wei);
    hpvd10->Fill(depolnos, wei);
    hpvd11->Fill(depolnos, wei);
    hpvd12->Fill(fracE, wei);
    hpvd13->Fill(depol, wei);
    hpvd14->Fill(depol, wei);
    hpvd15->Fill(depol-depolnos,wei);
  }

}

void hamcExptPVDIS::RunSummary() {

  //  cout << "Made it to the summary "<<endl;

  //  hamcExpt::RunSummary(0);

}
