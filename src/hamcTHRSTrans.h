#ifndef ROOT_hamcTHRSTrans
#define ROOT_hamcTHRSTrans

//  hamcTHRSTrans -- Wrapper for the THRSTrans library
//  which is the separately maintained software Seamus Riordan wrote.
//  R. Michaels, Dec 2017

#include "THRSTrans.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <map>

class hamcSpecHRS;
class THRSTrans;  // Seamus's class for 2nd order Transport. 

class hamcTHRSTrans : public hamcTrans {

  public:

     hamcTHRSTrans();
     virtual ~hamcTHRSTrans(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

  protected:

     THRSTrans *trans;
     THRSTrans::tune_t fTune; 

  private: 

     TString DoInterpret( THRSTrans::tune_t thetune);

     hamcTHRSTrans(const hamcTHRSTrans& trans);
     hamcTHRSTrans& operator=(const hamcTHRSTrans& trans);

#ifndef NODICT
ClassDef (hamcTHRSTrans, 0)   // Wrapper for the THRSTrans library
#endif
  };

#endif



   
