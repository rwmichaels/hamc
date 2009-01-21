//  hamcExpt   -- abstract base class for experiment
//  R. Michaels  Sept 2007

#include "hamcExpt.h"
#include "hamcPhysics.h"
#include "hamcSpecHRS.h"
#include "hamcTarget.h"
#include "hamcInout.h"
#include "hamcEvent.h"

using namespace std;

#ifndef NODICT
ClassImp(hamcExpt)
#endif



hamcExpt::hamcExpt(string sname): name(sname),
     numevent(0), numiter(1),
     setupfile("hamcsetup.dat"),
     didinit(kFALSE)
{
  physics = 0;    // must be implemented by inheriting class
  target  = 0;    //  "             "            "
  event   = 0;    //  "             "            " 
  inout   = new hamcInout();
  run_time = 720;  // hours
}

hamcExpt::~hamcExpt() {
  if (physics) delete physics;
  if (target)  delete target;
  if (event)   delete event;
  delete inout;
  for (std::vector<hamcSpecHRS*>::iterator ith = spectrom.begin();
      ith != spectrom.end(); ith++) delete *ith;
}

Int_t hamcExpt::InitInput(string sfile) {
  setupfile = sfile;
  inout->Init(this);
  return OK;
}

Int_t hamcExpt::Init(string sfile) {
  InitInput(sfile);
  return Init();
}

Int_t hamcExpt::Init() {

  inout->Init(this);  

  SetIterate(); 

  if(event) event->InitBeam(this);

  if(target) target->Init(this);

  if(physics) physics->Init(this);

  if(event) event->Init(this);

  inout->BookNtuple();

  hamcStrParser parser;
  parser.Load(inout->GetStrVect("run_time"));
  if (parser.IsFound("run_time")) {
     run_time = parser.GetData(); 
     cout << "Run time (hours) "<<run_time << endl;
  }   

  didinit = kTRUE;

  return OK;

}
 
Int_t hamcExpt::ChkSpIndex(Int_t idx) const {
  if (idx < 0 || idx > (Int_t)spectrom.size()-1) return ERROR;
  return OK;
}

Int_t hamcExpt::SetIterate() {
// Here we may do 2 iterations to obtain a 'kick' in something.
// The info about what to iterate is passed via setup file
// to the class that must kick.  
  numiter = inout->numiter;
  return OK;  
}


Int_t hamcExpt::Run(Int_t maxevent) {

  if (!didinit) {
     cout << "hamcExpt::Run:ERROR: Need to Init() first."<<endl;
     return ERROR;
  }

  for (iteration = 0; iteration < numiter; iteration++) {

    for (Int_t ievt = 0; ievt < maxevent; ievt++ ) {

      if (ievt > 0 && ((ievt%10000)==1)) cout << "event "<<ievt<<endl;
 
      if (event) event->Process(this);

    }
  }

  return OK;
  
}

void hamcExpt::RunSummary() {

  inout->Finish();

}

hamcSpecHRS* hamcExpt::GetSpectrom(Int_t i) const {
  if (i < 0 || i > GetNumSpectrom()) {
    // well, this is really bad
    cout << "hamcExpt::GetSpectrom:: FATAL INDEX ERROR"<<endl;
    cout << "Attempting to access spectrom # "<<i<<endl;
    cout << "Stopping ..."<<endl;
    exit(0);
  }
  return spectrom[i];
}


Int_t hamcExpt::ReadSetup() {

  // here want to do like THaInout::LoadFile
  // to make a vector of THaStrings.
  // Strings have keywords built in for target, spectrom, etc.

  // Determine number of iterations.

  for (Int_t i=0; i<(Int_t)setupstr.size(); i++) {
    //    if (string contains class:iterate:var, whiciter=string)
  }
 

  return OK;

}


