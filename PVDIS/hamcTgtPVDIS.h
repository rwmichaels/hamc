#ifndef ROOT_hamcTgtPVDIS
#define ROOT_hamcTgtPVDIS

//  hamcTgtPVDIS   -- the PVDIS LD2 target
//  D. Wang   R. Michaels   Jan 2009

#include "hamcTarget.h"
#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcTgtPVDIS : public hamcTarget {

// The PVDIS target

  public:

     hamcTgtPVDIS();
     virtual ~hamcTgtPVDIS();    
     Int_t Init(hamcExpt *exp);

  protected:


  private: 

     hamcTgtPVDIS& operator=(const hamcTgtPVDIS& tgt);
     hamcTgtPVDIS(const hamcTgtPVDIS& tgt);


#ifndef NODICT
ClassDef (hamcTgtPVDIS, 0)   // Target for PVDIS
#endif


};



#endif

   
