#ifndef ROOT_hamcTgtPREX
#define ROOT_hamcTgtPREX

//  hamcTgtPREX   -- the PREX Pb/C target 
//  R. Michaels  June 2008

#include "hamcTarget.h"
#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcTgtPREX : public hamcTarget {

// The PREX target

  public:

     hamcTgtPREX();
     virtual ~hamcTgtPREX();    
     Int_t Init(hamcExpt *exp);

  protected:


  private: 

     hamcTgtPREX& operator=(const hamcTgtPREX& tgt);
     hamcTgtPREX(const hamcTgtPREX& tgt);


#ifndef NODICT
ClassDef (hamcTgtPREX, 0)   // Target
#endif


};



#endif

   
