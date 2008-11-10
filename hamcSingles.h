#ifndef ROOT_hamcSingles
#define ROOT_hamcSingles

//  hamcSingles   -- base class for single-arm experiment
//  R. Michaels  June 2008

#include "hamcExpt.h"
#include <string>

class hamcSingles : public hamcExpt {

  public:

     hamcSingles(std::string sname);
     virtual ~hamcSingles();
     virtual Int_t Init(std::string sfile);
     virtual Int_t SetSpectrom(Int_t which, Float_t pmom, Float_t theta);
     void SetP0(Float_t p) { P0 = p; };
     void SetTheta(Float_t ang) { angle = ang; };

  protected:

     Float_t P0, angle;

  private: 

     hamcSingles(const hamcSingles& expt);
     hamcSingles& operator=(const hamcSingles& expt);

#ifndef NODICT
ClassDef (hamcSingles, 0)   // Single experiment
#endif

};

#endif



   
