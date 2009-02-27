//**********************************************************************
//
//   Main code for running hamc
//   The Hall A Monte Carlo
//   R. Michaels, version 0, Sept, 2008
//
//**********************************************************************

#include "hamcExptPREX.h"
#include "TROOT.h"
#include "TRint.h"
#include <signal.h>

using namespace std;

int main(int argc, char **argv) 
{

// We use hamcPREX for PREX.

  hamcExptPREX *prex;
  prex = new hamcExptPREX();

// Here's where you might instead do, perhaps depending on argv
//   dvcs = new hamcExptDVCS();
//      or
//   happex = new hamcExptHAPPEX();


  Int_t nevents = 200000;

  if (argc >= 2) {
    nevents = atoi(argv[1]);
  }
 
  cout << "Number of events to process "<<nevents<<endl;

  string setupfile="prex.dat";

  prex->Init(setupfile);

  prex->Run(nevents);

  return 1;

}


