#ifndef ROOT_hamcInout
#define ROOT_hamcInout

//  hamcInout   -- input/output for hamc
//  R. Michaels  May 2008

#include "Rtypes.h"
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include "TTree.h"
#include <fstream>


class hamcExpt;
class hamcEvent;
class THaString;
class TFile;
class TNtuple;

using namespace std;

class hamcHist {
// Utility class for histograms (see BookHisto method)
 public: 
  hamcHist(std::string nm, std::string title, Int_t bk, Float_t *dp,
           Int_t nb, Float_t xl, Float_t xh, Bool_t iterate, Bool_t acc) {
    name=nm; bkpoint=bk;  xdata = dp; ydata=0;
    dimen=1; weighted=kFALSE;
    care_accept=acc;
    fH2 = 0;  fH2x = 0;
    fH1 = new TH1F(name.c_str(), title.c_str(), nb, xl, xh);
    fH1x = 0;
    if (iterate) {
       name = name + "x";
       title = title + "_ITERATION";
       fH1x = new TH1F(name.c_str(), title.c_str(), nb, xl, xh);
    }
  };
  hamcHist(std::string nm, std::string title, Int_t bk, 
  	            Float_t *dpx,  Float_t *dpy,
                    Int_t nbx, Float_t xl, Float_t xh, 
                    Int_t nby, Float_t yl, Float_t yh, 
                    Bool_t iterate, Bool_t acc) {
    name=nm; bkpoint=bk; xdata = dpx; ydata = dpy;
    dimen=2; weighted=kFALSE;
    care_accept=acc;
    fH1 = 0;  fH1x = 0;
    fH2 = new TH2F(name.c_str(), title.c_str(), 
                 nbx, xl, xh, nby, yl, yh);
    fH2x = 0;
    if (iterate) {
      name = name + "x";
      title = title + "_ITERATION";
      fH2x = new TH2F(name.c_str(), title.c_str(), 
         nbx, xl, xh, nby, yl, yh);
    }
  };
  ~hamcHist() { 
     if (fH1) delete fH1;  
     if (fH1x) delete fH1x;  
     if (fH2) delete fH2;  
     if (fH2x) delete fH2x; 
  };
  void Fill(Int_t iter, Int_t bk, Float_t wei=1, Int_t inacc=1) {
// Fill histogram (1D or 2D / weighted or not / for 0th or 1st iteration) 
    if (bk != bkpoint) return;  // Only for this break point.
    if (care_accept && !inacc) return;  // If we care and if outside acceptance
    if (dimen==1) {    // 1D histogram
      if (weighted) {
         if (iter==0) {
           if (fH1) { fH1->Fill(*xdata,wei); } else { PrintErr(iter); }
	 } else {
           if (fH1x) { fH1x->Fill(*xdata,wei); } else { PrintErr(iter); }
	 }
      } else {
	if (iter==0) {
	   if (fH1) { fH1->Fill(*xdata); } else { PrintErr(iter); }
	} else {
	  if (fH1x) { fH1x->Fill(*xdata); } else { PrintErr(iter); }
	}
      }
    } else {     // 2D histogram
      if (weighted) {
        if (iter==0) {
          if (fH2) { fH2->Fill(*xdata,*ydata,wei); } else { PrintErr(iter); }
	} else {
          if (fH2x) { fH2x->Fill(*xdata,*ydata,wei); } else { PrintErr(iter); }
	}
      } else {
	if (iter==0) {
          if (fH2) { fH2->Fill(*xdata,*ydata); } else { PrintErr(iter); }
	} else {
          if (fH2x) { fH2x->Fill(*xdata,*ydata); } else { PrintErr(iter); }
	}
      }
    }
  };
  void SetWeighted() { weighted = kTRUE; };
  void Write() { 
    if (fH1)  fH1->Write();
    if (fH1x) fH1x->Write();
    if (fH2)  fH2->Write();
    if (fH2x) fH2x->Write();
  };
 private:  
  void PrintErr(Int_t iter) {
// This might mean you booked a histogram before finding out if there
// was an iteration, hence fH*x =0, and then you tried to do an iteration.
    std::cout << "hamcHist:ERROR: Inconsistent combination for iter "<<
       iter<<" "<<dimen<<std::endl;
  };
  TH1F *fH1, *fH1x;  // 1D, 1st and 2nd iteration(x)
  TH2F *fH2, *fH2x;  // 2D,   "   "
  std::string name;
  Int_t bkpoint, dimen;
  Float_t *xdata, *ydata;
  Bool_t weighted, care_accept;
};    

class hamcStrParser {
// Utility class to parse the strings, e.g. for iteration
//    "kick:track  theta  0.45"
 public:
  hamcStrParser() { idx = -1;  sdata.clear(); };
  void Load(std::vector<std::string> sd) { sdata = sd;  idx = -1; };
  void Print() {
    std::cout << "hamcStrParser:  "<<std::endl;
    std::cout << "idx = "<<idx<<"   sdata size "<<sdata.size()<<std::endl;
    if (sdata.size() > 0) {
      for (Int_t i=0; i<(Int_t)sdata.size(); i++) {
	std::cout << sdata[i]<<"   ";
      }
      std::cout<<std::endl;
    }
  }
  Bool_t IsFound(std::string str) {
    idx = -1;
    if ((Int_t)sdata.size() == 0) return kFALSE;
    for (Int_t i=0; i<((Int_t)sdata.size()); i++) {
      if (sdata[i] == str) {
   	idx = i+1;   // data is next word
        return kTRUE;
      }
    }
    return kFALSE;
  }
  Float_t GetData() {
    if (idx < 0 || idx > (Int_t)sdata.size()) return 0;
    Float_t result;
    sscanf(sdata[idx].c_str(),"%f",&result);
    return result;
  }
 private:
  std::vector<std::string> sdata;
  Int_t idx;
};
  
class hamcInout {

  public:

     hamcInout();
     virtual ~hamcInout(); 
 
     Int_t numiter;

     Int_t Init(hamcExpt*);  
     Int_t BookNtuple();
     Int_t Process(hamcExpt*);
     Int_t NumInMap(std::string);
     void SetWeight(Float_t wei) { weight = wei; };
     Float_t GetWeight() { return weight; };  
     Bool_t FoundString(std::string);
     std::vector<std::string> GetStrVect(std::string, Int_t i=0);
     Int_t AddToNtuple(std::string var, Float_t* dptr);
     Int_t AddToNtuple(std::string var, Int_t* dptr);
     TTree *t1;
     ofstream deriv;
     Int_t Finish();

     Int_t BookHisto(Bool_t tow, Bool_t acc, Int_t bk, std::string hid, std::string stitle, Float_t* dx, Int_t nxbin, Float_t xlo, Float_t xhi, Float_t* dy=0, Int_t nybin=0, Float_t ylo=0, Float_t yhi=0);

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

  protected:

     
  private: 

     Int_t LoadFile(std::string file);
     Int_t LoadFile();
     Int_t SetNumIterations();
     Int_t FillHisto(Int_t id, Float_t x, Float_t y=0, Float_t wei=1, Int_t iter=0);

     Bool_t did_init, root_disable, ntup_disable, isdone;
     Float_t weight;
     std::string setupfile;
     TFile *hFile;
     std::vector<hamcHist* > fHist;
     std::multimap<std::string, std::vector<std::string> > sdata;
     TNtuple *ntup;
     Float_t *fntup;
     std::string sntup;
     std::vector<Float_t *> fntptr;
     std::vector<Int_t *> intptr;

     hamcInout(const hamcInout& xout);
     hamcInout& operator=(const hamcInout& xout);


#ifndef NODICT
ClassDef (hamcInout, 0)   // input/output
#endif

};

#endif



   
