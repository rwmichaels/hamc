#ifndef ROOT_hamcExptPVDIS
#define ROOT_hamcExptPVDIS

//  hamcExptPVDIS   -- Pb Radius Experiment
//  D. Wang  R. Michaels  Jan. 2009

#include "hamcExpt.h"
#include "hamcSingles.h"
#include "TH1F.h"
#include "TH2F.h"

class hamcExptPVDIS : public hamcSingles {

  public:

     hamcExptPVDIS();
     virtual ~hamcExptPVDIS();
     Int_t Init(std::string sfile);

     void EventAnalysis();
     void RunSummary();

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptPVDIS(const hamcExptPVDIS& expt);
     hamcExptPVDIS& operator=(const hamcExptPVDIS& expt);

     TH1F *hpvd1, *hpvd2, *hpvd3;
     TH2F *hpvd4;
     TH1F *hpvd5, *hpvd6, *hpvd7;
     TH2F *hpvd8;
     TH1F *hpvd9, *hpvd10, *hpvd11;
     TH1F *hpvd12, *hpvd13, *hpvd14;
     TH1F *hpvd15;



#ifndef NODICT
ClassDef (hamcExptPVDIS, 0)   // PVDIS Experiment
#endif


};


#endif


   
