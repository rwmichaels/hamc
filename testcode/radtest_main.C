// Standalone code to check radcor code

#define NUMEVENT  1000000   // for testing

#include <iostream>
#include <string>
#include <iostream.h>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "hamcRad.h"

using namespace std;

int main(int argc, char* argv[])
{

// Radiation class
  hamcRad *rad = new hamcRad();

  Float_t E0 = 1.05;   // beam energy, GeV
  Float_t theta = 5;   // Central angle, degrees
  Float_t Zavg = 82;   // weighted Z (charge) of target
  Float_t trlen = 0.1;  // radiation length (fractional) of target
  Float_t tlen = 0.005;  // length (meters) of target

  rad->Init(E0, theta, Zavg, trlen, tlen);

// Initialize root and output
  TROOT fadcana("radroot","Radcor test");
  TFile hfile("radc.root","RECREATE","Radiative correction test");

  TH1F *hbint = new TH1F("hbint","Internal Brems.",1000,-0.1,1.04*E0);

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

  Float_t deI, deEin, deEout;

  for (Int_t iev = 0; iev < NUMEVENT; iev++) {

    Float_t tgt = tlen*gRandom->Rndm();

    rad->Generate(tgt);

    deI = rad->GetDeIntern();

    hbint->Fill(deI);

    Int_t idx = rad->LookupIdx(tgt);
 
    deEin = rad->GetDeExternIn();
    deEout = rad->GetDeExternOut();

    slice1[idx]->Fill(deEin);
    slice2[idx]->Fill(deEout);

  }

  hfile.Write();
  hfile.Close();

}

