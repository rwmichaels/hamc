#ifndef ROOT_hamcExptHAPPEX
#define ROOT_hamcExptHAPPEX

//  hamcExptHAPPEX   -- HAPPEX (eventually 3).
//  R. Michaels  Dec 2008

#include "hamcExpt.h"
#include "hamcSingles.h"

class hamcExptHAPPEX : public hamcSingles {

  public:

     hamcExptHAPPEX();
     virtual ~hamcExptHAPPEX();
     Int_t Init(std::string sfile);
     void Analysis();

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptHAPPEX(const hamcExptHAPPEX& expt);
     hamcExptHAPPEX& operator=(const hamcExptHAPPEX& expt);

#ifndef NODICT
ClassDef (hamcExptHAPPEX, 0)   // HAPPEX Experiment
#endif


};


#endif


   
