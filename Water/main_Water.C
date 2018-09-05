//**********************************************************************
//
//   Main code for running hamc
//   The Hall A Monte Carlo
//   R. Michaels, version 0, Sept, 2018
//
//**********************************************************************

#include "hamcExptWater.h"
#include "TROOT.h"
#include "TRint.h"
#include <signal.h>

using namespace std;

int main(int argc, char **argv) 
{


  hamcExptWater *watercell;
  watercell = new hamcExptWater();

// Here's where you might instead do, perhaps depending on argv
//   dvcs = new hamcExptDVCS();
//      or
//   watercell = new hamcExptWater();


  Int_t nevents = 200000;

  if (argc >= 2) {
    nevents = atoi(argv[1]);
  }
 
  cout << "Number of events to process "<<nevents<<endl;

  string setupfile="watercell.dat";

  watercell->Init(setupfile);

  watercell->Run(nevents);

  return 1;

}


