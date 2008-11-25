//  hamcPhysics   -- abstract base class for the physics
//  D. Jaunzeikare, R. Michaels  May 2008

#include "hamcPhysics.h"
#include "hamcExpt.h"
#include "hamcInout.h"
#include "hamcKine.h" 
#include "hamcRad.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifndef NODICT
ClassImp(hamcPhysics)
#endif


hamcPhysics::hamcPhysics(): did_init(kFALSE),crsec(0),asymmetry(0)
{
  radiation = new hamcRad();
  kine = new hamcKine();
}

hamcPhysics::~hamcPhysics()
{
  delete radiation;
  delete kine;
}

Int_t hamcPhysics::Init(hamcExpt* expt) {
// Here you want to grab from "expt" the parameters you need
// which depends on experiment and is the same for all events,
// i.e. definition of target

  expt->inout->AddToNtuple("crsec",&crsec);
  expt->inout->AddToNtuple("asy",&asymmetry);

  radiation->Init(expt);
  kine->Init(expt);

  did_init = kTRUE;
  return 1;
}


Int_t hamcPhysics::Radiate(hamcExpt* expt) {

  if (do_radiate) {
    return radiation->Generate(expt);
  } else {
    return 1;
  }
}

Int_t hamcPhysics::Generate(hamcExpt* expt) {

// This routine will be called for each event.
// First you'll need to extract from "expt" things which
// change each event, like energy and angles.

  CrossSection();

  Asymmetry();  

  return 1;
}


Int_t hamcPhysics::CrossSection() {
// Generates the cross section. Called by Generate.
  crsec = 0;
  return 1;
}

Int_t hamcPhysics::Asymmetry() {
// Generates the asymmetry. Called by Generate.
  asymmetry = 0;
  return 1;
}

