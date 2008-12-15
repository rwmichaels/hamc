#ifndef ROOT_hamcEloss
#define ROOT_hamcEloss

//  hamcEloss   -- 
//    Internal and external Brehmstrahlung.
//    Ionization loss.
//  R. Michaels  Nov 2008

#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcEloss {

  public:

     hamcEloss();
     virtual ~hamcEloss(); 

     Int_t Init(hamcExpt *expt);
     Int_t InitRad(Float_t E,Float_t theta,Float_t z,Float_t rl,Float_t tl);

     Int_t Generate(hamcExpt *expt);  // event generator
     Int_t GenerateRad(Float_t rlin, Float_t rlout); 
     Int_t GenerateRad(Float_t zpos);
     Int_t GenerateDeDx(hamcExpt *expt);
     Float_t GetDeIntern();  
     Float_t GetDeExternIn();  
     Float_t GetDeExternOut();  
     Float_t GetDeIonIn() { return dE_IonizeIn; };
     Float_t GetDeIonOut() { return dE_IonizeOut; };
// dE/dx from collisional ionization
     Float_t dedx_eloss(Float_t znuc);
     Int_t LookupIdx(Float_t tl);
     Int_t CheckInit();
     void Print();
// Methods taken from gener_cone, for tests (option for use)
     Float_t gener_radlossint(Float_t k, Float_t nu);
     Float_t gener_radlossext(Float_t k, Float_t dist_tgt);

  private:

     void Setup(Int_t which, Float_t tl);

     static const Int_t MAXCNT=300000;
     static const Int_t ldebug=0;
     static const Float_t Me=0.0000511;
     static const Float_t Euler=0.5772157;
     static const Float_t fracresol=1e-4;

     Bool_t did_init;
     Bool_t use_genercone, use_ionize;
     Int_t Npts,nybin,Nslices;
     Float_t tequiv,trlen,tlen,tgtZ,tdensity,E0,qsq,dE;
     Float_t me,alpha,pi,yfact,ycell,bval;
     Float_t dE_IntBrehm, dE_ExtBrehmIn, dE_ExtBrehmOut, dE_Ionization;
     Float_t dE_Bsum;
     Float_t dE_IonizeIn, dE_IonizeOut; 
     std::vector<std::vector<Float_t> > Eextern;
     std::vector<Float_t > Eintern;
// Energy resolution cuts
     Float_t cut_intern;
     std::vector<Float_t > cut_extern;

     hamcEloss(const hamcEloss& eloss);
     hamcEloss& operator=(const hamcEloss& eloss);


#ifndef NODICT
ClassDef (hamcEloss, 0)   // electromagnetic radiation
#endif

};

#endif




   
