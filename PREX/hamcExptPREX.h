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
     void RunSummary(Int_t i=0);

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptPREX(const hamcExptPREX& expt);
     hamcExptPREX& operator=(const hamcExptPREX& expt);

     TH2F *prex_xy1, *prex_xy2, *prex_xy3, *prex_xy4, *prex_xy5;
     TH1F *prex_x1, *prex_x2, *prex_x2a, *prex_x3;

     Float_t sumr_pc1, xcnt_pc1, sumr_pc2, xcnt_pc2;
     Float_t solid_athole;
     Float_t inatdet;
     Float_t xdetlo, xdethi, ydetlo, ydethi; // A_T detector


#ifndef NODICT
ClassDef (hamcExptPREX, 0)   // PREX Experiment
#endif


};


#endif


   
