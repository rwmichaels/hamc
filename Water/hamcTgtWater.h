#ifndef ROOT_hamcTgtWater
#define ROOT_hamcTgtWater

//  hamcTgtWater   -- the Water Cell target
//  R. Michaels  Dec 2018

#include "hamcTarget.h"
#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcTgtWater : public hamcTarget {

// The Water cell target

  public:

     hamcTgtWater();
     virtual ~hamcTgtWater();    
     Int_t Init(hamcExpt *exp);

  protected:


  private: 

     hamcTgtWater& operator=(const hamcTgtWater& tgt);
     hamcTgtWater(const hamcTgtWater& tgt);


#ifndef NODICT
ClassDef (hamcTgtWater, 0)   // Target for Water Cell
#endif


};



#endif

   
