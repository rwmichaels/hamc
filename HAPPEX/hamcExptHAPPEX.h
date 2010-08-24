#ifndef ROOT_hamcExptHAPPEX
#define ROOT_hamcExptHAPPEX

//  hamcExptHAPPEX   -- HAPPEX (eventually 3).
//  R. Michaels  Dec 2008

#include "hamcExpt.h"
#include "hamcSingles.h"
#include "TH1F.h"
#include "TH2F.h"

class hamcExptHAPPEX : public hamcSingles {

  public:

     hamcExptHAPPEX();
     virtual ~hamcExptHAPPEX();
     Int_t Init(std::string sfile);
     void EventAnalysis();

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptHAPPEX(const hamcExptHAPPEX& expt);
     hamcExptHAPPEX& operator=(const hamcExptHAPPEX& expt);

     TH1F *qsq1, *qsq2, *qsq3, *qsq4;
     TH2F *hxy1, *hxy2, *hxy3, *hxy4;


#ifndef NODICT
ClassDef (hamcExptHAPPEX, 0)   // HAPPEX Experiment
#endif


};


#endif


   
