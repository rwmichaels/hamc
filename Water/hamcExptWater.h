#ifndef ROOT_hamcExptWater
#define ROOT_hamcExptWater

//  hamcExptWater   -- Water Cell
//  R. Michaels  Aug 2018

#include "hamcExpt.h"
#include "hamcSingles.h"
#include "TH1F.h"
#include "TH2F.h"

class hamcExptWater : public hamcSingles {

  public:

     hamcExptWater();
     virtual ~hamcExptWater();
     Int_t Init(std::string sfile);
     void EventAnalysis();
     Float_t O16CrossSection(Float_t energy, Float_t angle);

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptWater(const hamcExptWater& expt);
     hamcExptWater& operator=(const hamcExptWater& expt);

     TH1F *qsq1, *qsq2, *qsq3, *qsq4;
     TH2F *hxy1, *hxy2, *hxy3, *hxy4;


#ifndef NODICT
ClassDef (hamcExptWater, 0)   // Water Cell Experiment
#endif


};


#endif


   
