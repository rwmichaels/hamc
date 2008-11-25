// Standalone code to checkhamcKine

#define NUMEVENT  100000   // for testing

#include <iostream>
#include <string>
#include <iostream.h>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "hamcKine.h"
#include "hamcRad.h"

using namespace std;

int main(int argc, char* argv[])
{

// Kinematics class
  hamcKine *kine = new hamcKine();

// Radiation class
  hamcRad *rad = new hamcRad();

  Float_t E0 = 1.05;    // beam energy, GeV
  Float_t theta = 5;   // Central angle, degrees
  Float_t pi = 3.1415926;
  Float_t thetamin = 2*pi/180;
  Float_t thetamax = 8*pi/180;
  Float_t Zavg = 82;    // weighted Z (charge) of target
  Float_t mass_tgt = 208*0.938;  // mass of target (GeV)
  Float_t trlen = 0.1;  // radiation length (fractional) of target
  Float_t tlen = 0.005;  // length (meters) of target  
  Float_t phimin = 0;
  Float_t phimax = pi;
 
  string process = "elastic";    // "elastic" or "dis" supported

  if (process == "elastic") {
    kine->Init(process, E0, theta, mass_tgt, 
	     thetamin, thetamax, phimin, phimax);   
  } else {
    kine->Init(process, E0, theta, mass_tgt, 
	     thetamin, thetamax, phimin, phimax, 0.1*E0, 0.9*E0);   
    kine->SetDisDef(0.1,0.99,2,4);
  }

  rad->Init(E0, theta, Zavg, trlen, tlen);

// Initialize root and output
  TROOT hamcana("kinroot","Kinematics generator test");
  TFile hfile("kine.root","RECREATE","Kinematics generation test");

  TH1F *he1 = new TH1F("he1","Energy generated",1000,0.2*E0,1.2*E0);
  TH1F *he2 = new TH1F("he2","Eprime",100,0.2*E0,1.2*E0);
  TH2F *hepe = new TH2F("hepe","Eprime vs energy ",
 			 100,0.2*E0,1.2*E0,
			 100,0.2*E0,1.2*E0);
  TH2F *hepth = new TH2F("hepth","Eprime vs energy ",100,
                         0.8*thetamin,1.2*thetamax,
			 100,0.2*E0,1.2*E0);
  TH2F *hthph = new TH2F("hthph","Phi vs Theta",100,
             0.8*thetamin,1.2*thetamax,100,0.8*phimin-pi/8,1.2*phimax);

  TH1F *hqsq;  TH1F *hmass; TH2F *hwvq;
  if (process == "elastic") {
    hqsq = new TH1F("hqsq","Qsq distribution",100,-0.02,0.15);
    hmass = new TH1F("hmass","Inv. Mass ",100,-2,50);
    hwvq = new TH2F("hwvq","Wsq vs qsq ",100,-0.1,0.15,100,-2,50);
  } else {
    hqsq = new TH1F("hqsq","Qsq distribution",100,0.1,15);
    hmass = new TH1F("hmass","Inv. Mass ",100,-0.1,5);
    hwvq = new TH2F("hwvq","Wsq vs qsq ",100,1,10,100,2,15);
  }
  

  TH1F *hbint = new TH1F("hbint","Internal Brems.",1000,-0.1,1.04*E0);
  TH1F *hbext = new TH1F("hbext","External Brems.",1000,-0.1,1.04*E0);

  vector<TH1F *> slice1, slice2;
  char cname[50],ctitle[100];

  for (Int_t i=0; i<10; i++) {
    sprintf(cname,"h1slice%d",i+1);
    sprintf(ctitle,"Straggling(in) for slice %d",i+1);
    slice1.push_back(new TH1F(cname,ctitle,1000,-0.1,1.04*E0));
    sprintf(cname,"h2slice%d",i+1);
    sprintf(ctitle,"Straggling(out) for slice %d",i+1);
    slice2.push_back(new TH1F(cname,ctitle,1000,-0.1,1.04*E0));
 }

  Float_t deI, dEin, dEout;
  Float_t ebeam;

  for (Int_t iev = 0; iev < NUMEVENT; iev++) {

    Float_t tgt = tlen*gRandom->Rndm();

    rad->Generate(tgt);

    deI = rad->GetDeIntern();

    hbint->Fill(deI);

    Int_t idx = rad->LookupIdx(tgt);
 
    dEin = rad->GetDeExternIn();
    dEout = rad->GetDeExternOut();

    hbext->Fill(dEin+dEout);

    slice1[idx]->Fill(dEin);
    slice2[idx]->Fill(dEout);

    ebeam = E0 * (1 + 1e-3*gRandom->Rndm());
    ebeam = ebeam - dEin;

    Int_t status = kine->Generate(ebeam, deI+dEout);
    if (status == -1) continue;

    he1->Fill(kine->energy);
    hthph->Fill(kine->theta,kine->phi); 
    hepe->Fill(kine->energy,kine->eprime);
    hepth->Fill(kine->theta,kine->eprime);
    he2->Fill(kine->eprime);
    hqsq->Fill(kine->qsq);
    hwvq->Fill(kine->qsq,kine->wsq);
    
    Float_t wsq = kine->wsq;
    Float_t xmass = 0;
    if (wsq >= 0) xmass = TMath::Sqrt(wsq);
    
    hmass->Fill(xmass);

  }

  hfile.Write();
  hfile.Close();

}

