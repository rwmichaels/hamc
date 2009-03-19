//  hamcInout   -- Input/Output including
//     1) Setup file (ASCII)
//     2) Histograms, written to a ROOT file
//     3) Ntuple, written to a ROOT file
//
//  R. Michaels  May 2008

#include "hamcInout.h"
#include "hamcExpt.h"
#include "hamcEvent.h"
#include "THaString.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TRegexp.h"
#include "TROOT.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcInout)
#endif


hamcInout::hamcInout() : numiter(0),did_init(kFALSE),weight(1),setupfile("hamc.dat"),hFile(0),fntup(0) {
  ntup  = 0;
  sntup = "";
}

hamcInout::~hamcInout() {
  if (fntup) delete [] fntup;
  if (hFile) delete hFile;
  for (std::vector<hamcHist*>::iterator ith = fHist.begin();
      ith != fHist.end(); ith++) delete *ith;
}


Int_t hamcInout::Init(hamcExpt *expt) {

// Initialize ROOT and the input/output.
// Call this *BEFORE* initializing other classes
// thta need I/O.

   if (did_init) return OK;  // already init'd 

   TROOT fadcana("hamcroot","Hall A Monte Carlo");
   hFile = new TFile("hamc.root","RECREATE","Hall A Monte Carlo");

   Int_t errcond = 0;
   setupfile = expt->GetSetupFile();
   cout << "Inout: setupfile "<<setupfile<<endl;

   if (LoadFile()==ERROR) {
     cout << "hamcInout:Init:ERROR"<<endl;
     errcond = 1;
   }
 
   SetNumIterations();

   did_init = kTRUE;

   return OK;

}

Int_t hamcInout::SetNumIterations() {
   numiter = 1;
   if (FoundString("iterate")) numiter = 2;
   return OK;
}


Int_t hamcInout::Process(hamcExpt *expt) {

// Fill histograms and ntuples .
// Must be called for each break point in spectrometer.
//
// Histogram filling depends on break point and iteration
// and whether you're in acceptance (and whether you care), see inacc.
// ntuple filled for 0th iteration only and only at focal plane.

  Int_t iter, ibrk, inacc;

  iter  = expt->iteration;        // iteration of the experiment
  ibrk  = expt->event->brkpoint;  // break point in the event loop
  inacc = expt->event->inaccept;  // in acceptance (1) or not (0)

// Flag "inacc" is 1 if we're in the acceptance, 0 if not.
// This may or may not affect histograms, depending on if
// they were booked to "care" about acceptance.

  for (Int_t id=0; id<(Int_t)fHist.size(); id++) {

    fHist[id]->Fill(iter, ibrk, weight, inacc);

  }


// Ntuple filled at focal point, for iteration=0 only.

  if (ibrk == IFOCAL) {

    if (iter == 0 && fntup !=0  && ntup !=0 ) {


      for (Int_t ivar = 0; ivar < (Int_t)fntptr.size(); ivar++) {

        if (fntptr[ivar] != 0) {     // floating
          fntup[ivar] = *fntptr[ivar];
	} else {                     // integer (same size vector)
          fntup[ivar] = *intptr[ivar];
	}

      }    
      ntup->Fill(fntup);

    }
  }
  
  return OK;
 
}

Int_t hamcInout::Finish() {
  
  if (hFile) {
    hFile->cd();
    for (Int_t i=0; i<(Int_t)fHist.size(); i++) {
       fHist[i]->Write();
    }
    if (ntup) ntup->Write();
    hFile->Write();
    hFile->Close();
  }
  return OK;
}

Int_t hamcInout::LoadFile(string sfile) {
  setupfile=sfile;
  return LoadFile();

}

Int_t hamcInout::LoadFile() {
  ifstream* odef = new ifstream(setupfile.c_str());
  if ( ! (*odef) ) {
    cout << "ERROR: hamcInout: Loadfile:  file "<<setupfile<<" not found !"<<endl;
    return ERROR;
  }
  const string comment = "#";
  const string whitespace( " \t" );
  string::size_type pos;
  vector<THaString> strvect;
  THaString sline;
  while (getline(*odef,sline)) {
    // Blank line or comment line?
    if( sline.empty()
	|| (pos = sline.find_first_not_of( whitespace )) == string::npos
	|| comment.find(sline[pos]) != string::npos )
      continue;
    // Get rid of trailing comments
    if( (pos = sline.find_first_of( comment )) != string::npos )
      sline.erase(pos);
    // Split the line into tokens separated by whitespace
    strvect = sline.Split();
    vector<string> srest;
    srest.clear();    
    if (strvect.size() > 1) {
      for (Int_t i=1; i<(Int_t)strvect.size(); i++) {
        srest.push_back(strvect[i]);
      }
    }
    sdata.insert(make_pair(strvect[0],srest));
  }

  return OK;
}


Int_t hamcInout::NumInMap(string key) {
  // return size of data in sdata for that key.
   return sdata.count(key);
}

Bool_t hamcInout::FoundString(string key) {
  return (NumInMap(key) > 0);
}

vector<string> hamcInout::GetStrVect(string key, Int_t i) {
// Get element #i from the sdata map.
// Although this is rather inefficient & slow, it's ok since 
// it's only done in the initialization phase.

  Int_t ncnt=0;
  vector<string> result;
  result.clear();

  multimap<string, vector<string> >::const_iterator lb = sdata.lower_bound(key);
  multimap<string, vector<string> >::const_iterator ub = sdata.upper_bound(key);
  for ( multimap<string, vector<string> >::const_iterator pm = lb; pm != ub; pm++) {
    if (i == ncnt++) return pm->second;
  }

  return result;

}

Int_t hamcInout::BookHisto(Bool_t toweight, Bool_t accept, Int_t bkpt, string hid, string stitle, Float_t *ptrx, Int_t nxbin, Float_t xlo, Float_t xhi, Float_t *ptry, Int_t nybin, Float_t ylo, Float_t yhi) {

// Arguments
//     toweight = logical, to weight or not (a typ. weight = cross section)
//     accept = logical, to use acceptance cut or not (if not, histo filled anyway)
//     bkpt = break point, where to fill histogram (e.g. ITARGET, ICOLLIM, etc)
//                 these are defined in hamcSpecHRS.h
//     hid = e.g. h5 and the histogram appears as "h5" in output (h5->Draw, etc)
//                 hid is a unique ID.
//     stitle = string title of histogram
//     ptrx = pointer to data to put in histogram.  This is a pointer to the
//              raw data in the class that calls BookHisto.
//     nxbin = num bins in X
//     xlo, xhi = range in X
//     similarly y (if nybin=0, it's a 1D histo)
//     More histos like "hid_ITERATE" appear if do_iterate is true.

  Int_t dimen = 1;
  if (nybin > 0) dimen=2;

  Bool_t do_iterate = kFALSE;
  if (numiter > 1) do_iterate = kTRUE;

  if (dimen==1) {
    fHist.push_back(new hamcHist(hid,stitle,bkpt,ptrx,nxbin,xlo,xhi,do_iterate,accept));
  } else {
    fHist.push_back(new hamcHist(hid,stitle,bkpt,ptrx,ptry,nxbin,xlo,xhi,nybin,ylo,yhi,do_iterate,accept));
  }

  Int_t idx = fHist.size()-1;
  if (toweight && idx >= 0) {
      fHist[idx]->SetWeighted();
  }

  return OK;

}

Int_t hamcInout::AddToNtuple(std::string var, Int_t* dptr) {
// Add integer data to ntuple

  sntup += var;
  sntup += ":";
  
  fntptr.push_back(0);  // keep indices in synch
  intptr.push_back(dptr);

  return OK;


}

Int_t hamcInout::AddToNtuple(std::string var, Float_t* dptr) {
// Add floating point data to ntuple

  sntup += var;
  sntup += ":";
  
  fntptr.push_back(dptr);
  intptr.push_back(0);  // keep indices in synch

  return OK;

}

Int_t hamcInout::BookNtuple() {

// Call this after all calls to "AddToNtuple", normally done in hamcEvent

  if (fntptr.size() == 0) return OK;

  string::size_type pos = sntup.find_last_of(":");
  string stemp = sntup.substr(0,pos);

  ntup  = new TNtuple("hamc","hamc ntuple",stemp.c_str());
  fntup = new Float_t[fntptr.size()+1];    
  return OK;

}
   



