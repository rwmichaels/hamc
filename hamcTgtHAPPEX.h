#ifndef ROOT_hamcTgtHAPPEX
#define ROOT_hamcTgtHAPPEX

//  hamcTgtHAPPEX   -- the HAPPEX LH2 target
//  R. Michaels  Dec 2008

#include "hamcTarget.h"
#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcTgtHAPPEX : public hamcTarget {

// The HAPPEX target

  public:

     hamcTgtHAPPEX();
     virtual ~hamcTgtHAPPEX();    
     Int_t Init(hamcExpt *exp);

  protected:


  private: 

     hamcTgtHAPPEX& operator=(const hamcTgtHAPPEX& tgt);
     hamcTgtHAPPEX(const hamcTgtHAPPEX& tgt);


#ifndef NODICT
ClassDef (hamcTgtHAPPEX, 0)   // Target for HAPPEX
#endif


};



#endif

   
