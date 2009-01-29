//  hamcPhyPVDIS   -- class for the physics of PVDIS
//  D. Jaunzeikare, R. Michaels  May 2008


#include "hamcPhyPVDIS.h"
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
ClassImp(hamcPhyPVDIS)
#endif


hamcPhyPVDIS::hamcPhyPVDIS() : hamcPhysics()
{
  phy_name = "PVDIS physics";
  scatt_process = "dis";
  do_radiate = kTRUE;
}


hamcPhyPVDIS::~hamcPhyPVDIS() { }


Int_t hamcPhyPVDIS::Init(hamcExpt* expt) {

  didinit = kTRUE;

  hamcPhysics::Init(expt);


  return 1;
}


Int_t hamcPhyPVDIS::Generate(hamcExpt *expt) {

   Double_t energy = expt->physics->kine->energy;
   Double_t theta = expt->physics->kine->theta;
   Double_t eprime = expt->physics->kine->eprime;
   Double_t xbj = expt->physics->kine->x;
   Double_t qsq = expt->physics->kine->qsq;

   Double_t asymm=0.0;

   Z=1.0;
   A=2.0;

   //   cout<<"Z="<<Z<<"\tA="<<A<<"\tE="<<energy<<"\ttheta="<<theta<<"\teprime="<<eprime<<endl;

// Compute the cross section for PVDIS
   crsec =  cross_section__(&Z,&A,&energy,&theta,&eprime);

// Compute the PV asymmetry 
   getpdf_mrst2003c__(&energy, &xbj, &qsq, &asymm);
   asymmetry = float(asymm);
   if(crsec!=crsec) crsec=0.0;

   cout<<"energy="<<energy<<"\teprime="<<eprime<<"\txsection="<<crsec<<"\tasymm="<<asymmetry<<endl;

   return OK;
}


