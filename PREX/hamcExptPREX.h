#ifndef ROOT_hamcExptPREX
#define ROOT_hamcExptPREX

//  hamcExptPREX   -- Pb Radius Experiment
//  R. Michaels  May 2008

#include "TH1F.h"
#include "TH2F.h"

#include "hamcExpt.h"
#include "hamcSingles.h"

class hamcExptPREX : public hamcSingles {

  public:

     hamcExptPREX();
     virtual ~hamcExptPREX();
     Int_t Init(std::string sfile);

     void EventAnalysis();
     void RunSummary();

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptPREX(const hamcExptPREX& expt);
     hamcExptPREX& operator=(const hamcExptPREX& expt);

     TH2F *prex_xy1, *prex_xy2, *prex_xy3;
     TH1F *prex_x1, *prex_x2;

     Float_t sumr_pc, xcnt_pc;


#ifndef NODICT
ClassDef (hamcExptPREX, 0)   // PREX Experiment
#endif


};


#endif


   
