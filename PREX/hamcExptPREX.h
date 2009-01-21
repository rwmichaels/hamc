#ifndef ROOT_hamcExptPREX
#define ROOT_hamcExptPREX

//  hamcExptPREX   -- Pb Radius Experiment
//  R. Michaels  May 2008

#include "hamcExpt.h"
#include "hamcSingles.h"

class hamcExptPREX : public hamcSingles {

  public:

     hamcExptPREX();
     virtual ~hamcExptPREX();
     Int_t Init(std::string sfile);

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptPREX(const hamcExptPREX& expt);
     hamcExptPREX& operator=(const hamcExptPREX& expt);

#ifndef NODICT
ClassDef (hamcExptPREX, 0)   // PREX Experiment
#endif


};


#endif


   
