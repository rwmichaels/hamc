#ifndef ROOT_hamcTrans
#define ROOT_hamcTrans

//  hamcTrans   -- abstract base class for transport
//  R. Michaels  June 2008

#include "Rtypes.h"
#include <string>

class hamcTransVect;
class hamcSpecHRS;
class hamcTrack;

class hamcTrans {

  public:

     hamcTrans();
     virtual ~hamcTrans()=0; 
     virtual void Print()=0;
     virtual Int_t Init(hamcSpecHRS *spec)=0;
// Transform a track from its origin to "where"
     virtual Int_t TransForm(hamcTrack *trk, Int_t where) const=0;  
// Go a drift distance 'dist', modifies tvec
     Int_t Drift(Float_t dist, hamcTrack *trk) const;

  protected:

     Bool_t did_init;
     Int_t which_spectrom;
     Float_t collim_distance;

  private: 

     hamcTrans(const hamcTrans& trans);
     hamcTrans& operator=(const hamcTrans& trans);

#ifndef NODICT
ClassDef (hamcTrans, 0)   // Transport
#endif

};

#endif



   
