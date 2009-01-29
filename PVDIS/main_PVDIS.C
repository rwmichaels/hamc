//**********************************************************************
//
//   Main code for running hamc
//   The Hall A Monte Carlo
//   R. Michaels, version 0, Sept, 2008
//
//**********************************************************************

#include "hamcExptPVDIS.h"
#include "TROOT.h"
#include "TRint.h"
#include <signal.h>

using namespace std;

int main(int argc, char **argv) 
{


  hamcExptPVDIS *pvdis;
  pvdis = new hamcExptPVDIS();

// Here's where you might instead do, perhaps depending on argv
//   dvcs = new hamcExptDVCS();
//      or
//   pvdis = new hamcExptPVDIS();


  Int_t nevents = 200000;

  if (argc >= 2) {
    nevents = atoi(argv[1]);
  }
 
  cout << "Number of events to process "<<nevents<<endl;

  string setupfile="pvdis.dat";

  pvdis->Init(setupfile);

  pvdis->Run(nevents);

  pvdis->RunSummary();

  return 1;

}


