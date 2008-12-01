#ifndef ROOT_hamcRad
#define ROOT_hamcRad

//  hamcRad   -- Internal and external Brehmstrahlung.
//  Eventually want to add ionization loss.
//  R. Michaels  Nov 2008

#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcRad {

  public:

     hamcRad();
     virtual ~hamcRad(); 

     Int_t Init(hamcExpt *expt);
     Int_t Init(Float_t E,Float_t theta,Float_t z,Float_t rl,Float_t tl);

     Int_t Generate(hamcExpt *expt);  // event generator
     Int_t Generate(Float_t ztgt); 
     Float_t GetDeIntern();  
     Float_t GetDeExternIn();  
     Float_t GetDeExternOut();  
     Int_t LookupIdx(Float_t tl);
     Int_t CheckInit();
     void Print();
 // Stuff taken from gener_cone, for tests 
     Float_t gener_radlossint(Float_t k, Float_t nu);
     Float_t gener_radlossext(Float_t k, Float_t dist_tgt);

  private:

     void Setup(Int_t which, Float_t tl);

     static const Int_t MAXCNT=100000;
     static const Int_t ldebug=0;
     static const Float_t Me=0.0000511;
     static const Float_t Euler=0.5772157;

     Bool_t did_init;
     Bool_t use_genercone;
     Int_t Npts,nybin;
     Int_t Nslices;
     Float_t tequiv,trlen,tlen,E0,qsq,dE;
     Float_t me,alpha,pi,yfact,ycell,bval;
     Float_t dE_IntBrehm, dE_ExtBrehmIn, dE_ExtBrehmOut, dE_Ionization;
     Float_t dE_Bsum;
     std::vector<std::vector<Float_t> > Estraggle;
     std::vector<Float_t > Eintern;

     hamcRad(const hamcRad& phys);
     hamcRad& operator=(const hamcRad& phys);


#ifndef NODICT
ClassDef (hamcRad, 0)   // electromagnetic radiation
#endif

};

#endif




   
