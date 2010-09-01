#ifndef ROOT_hamcSingles
#define ROOT_hamcSingles

//  hamcSingles   -- base class for single-arm experiment
//  R. Michaels  June 2008

#include "hamcExpt.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>
#include <vector>

class hamcAccAvg;

class hamcSingles : public hamcExpt {

  public:

     hamcSingles(std::string sname);
     virtual ~hamcSingles();
     virtual Int_t Init(std::string sfile);
     virtual Int_t SetSpectrom(Int_t which, Float_t pmom, Float_t theta);
     void SetP0(Float_t p) { P0 = p; };
     void SetTheta(Float_t ang) { angle = ang; };
     virtual void EventAnalysis();
     virtual void RunSummary(Int_t iter=0);
     Int_t Run(Int_t maxevent);

  protected:

     Float_t P0, angle;
     std::vector<hamcAccAvg* > acc;
     Int_t num_mtl, num_phyt;
     Float_t dpp_cut;

  private: 

     hamcSingles(const hamcSingles& expt);
     hamcSingles& operator=(const hamcSingles& expt);

#ifndef NODICT
ClassDef (hamcSingles, 0)   // Single experiment
#endif

};

#endif



   
