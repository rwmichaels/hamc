#ifndef ROOT_hamcBeam
#define ROOT_hamcBeam

//  hamcBeam   -- the beam track
//  R. Michaels  Apr 2008

#include "Rtypes.h"
#include "hamcTrack.h"
#include "TLorentzVector.h"
#include <vector>
#include <map>

class hamcExpt;

class hamcBeam : public hamcTrack {

  public:

     hamcBeam();
     hamcBeam(Float_t ebeam, Float_t esigma);
     virtual ~hamcBeam(); 
   
     Int_t Init(hamcExpt *exp);
     Int_t Generate(hamcExpt *exp);
     Bool_t IsRastered() const { return rastered; };
     Int_t SetRaster(Float_t xmax, Float_t ymax);
     Float_t GetRaster(Int_t ixy);
     Int_t Radiate(hamcExpt *exp);
     void Print();

  protected:

  private: 

     Bool_t rastered;
     Float_t xrast,yrast;

     hamcBeam(const hamcBeam& phys);
     hamcBeam& operator=(const hamcBeam& phys);

#ifndef NODICT
ClassDef (hamcBeam, 0)   // Incoming beam
#endif


};

#endif



   
