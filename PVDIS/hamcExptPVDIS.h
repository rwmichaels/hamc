#ifndef ROOT_hamcExptPVDIS
#define ROOT_hamcExptPVDIS

//  hamcExptPVDIS   -- Pb Radius Experiment
//  D. Wang  R. Michaels  Jan. 2009

#include "hamcExpt.h"
#include "hamcSingles.h"

class hamcExptPVDIS : public hamcSingles {

  public:

     hamcExptPVDIS();
     virtual ~hamcExptPVDIS();
     Int_t Init(std::string sfile);

  private: 

  // Copy constructor and operator= defined null and private
     hamcExptPVDIS(const hamcExptPVDIS& expt);
     hamcExptPVDIS& operator=(const hamcExptPVDIS& expt);

#ifndef NODICT
ClassDef (hamcExptPVDIS, 0)   // PVDIS Experiment
#endif


};


#endif


   
