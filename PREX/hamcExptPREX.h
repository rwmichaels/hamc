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

     TH1F *prex_th0, *prex_th1;
     TH2F *prex_xy0, *prex_xy1, *prex_xy2, *prex_xy3, *prex_xy4, *prex_xy5;
     TH1F *prex_x1, *prex_x2, *prex_x2a, *prex_x3;
     TH1F *prex_theta, *prex_indet;
     TH1F *hqsq1, *hqsq2, *hqsq3, *hqsq4, *hqsq5, *hqsq6, *hqsq7, *hqsq8;
     TH1F *hqsq00, *hqsq10,*hqsq20,*hqsq30,*hqsq40,*hqsq50,*hqsq60,*hqsq70;
     TH1F *hqsq01,*hqsq02,*hqsq03,*hqsq04,*hqsq05,*hqsq06,*hqsq07;
     TH1F *hph00, *hph10,*hph20,*hph30,*hph40,*hph50,*hph60,*hph70;
     TH1F *hph01,*hph02,*hph03,*hph04,*hph05,*hph06,*hph07;
     TH1F *hth00, *hth10,*hth20,*hth30,*hth40,*hth50,*hth60,*hth70;
     TH1F *hth01,*hth02,*hth03,*hth04,*hth05,*hth06,*hth07;
     TH1F *hqsqmid;

     TH1F *qsqf;

     TH1F *hpx1;

     Float_t sumr_pc1, xcnt_pc1, sumr_pc2, xcnt_pc2;
     Float_t solid_athole;
     Float_t inatdet;
     Float_t xdetlo, xdethi, ydetlo, ydethi; // A_T detector


#ifndef NODICT
ClassDef (hamcExptPREX, 0)   // PREX Experiment
#endif


};


#endif


   
