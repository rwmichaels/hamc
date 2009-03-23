#ifndef ROOT_hamcEloss
#define ROOT_hamcEloss

//  hamcEloss   -- 
//    Internal and external Brehmstrahlung.
//    Ionization loss.
//  R. Michaels  Nov 2008

#include "Rtypes.h"
#include "TF1.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcEloss {

  public:

     hamcEloss();
     virtual ~hamcEloss(); 

     Int_t Init(hamcExpt *expt);
     Int_t InitRad();

     Int_t Generate(hamcExpt *expt);  // event generator
     Int_t Generate_gcone();
     Int_t GenerateNumer();
     Int_t Generate_tf1();
     Int_t GenerateDeDx(hamcExpt *expt);
     Float_t GetDeDx(Float_t radlen, Int_t where);

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

     void Setup_numer(Int_t which, Float_t tl);
     void Setup_tf1(Int_t which, Double_t rad, Double_t Ztgt, Double_t Ene); 
     Double_t GenLandLo(Float_t lambda, Float_t radlen);
     Double_t GenLandHi(Float_t lambda, Float_t radlen);

     static const Int_t MAXCNT=500000;
     static const Int_t ldebug=0;
     static const Float_t Me=0.0000511;
     static const Float_t Euler=0.5772157;
     static const Double_t pi=3.1415926;

     Bool_t did_init;
     Bool_t use_genercone, use_tf1, use_numer;
     Bool_t use_ionize;
     Int_t Npts,nybin,Nslices;
     Float_t psi_scale; 
     Float_t trlen,tequiv;  // RL and equiv. raditor (int. Brehm)
     Float_t tlen,tgtA,tgtZ,tdensity,E0,theta_central,qsq;
     Float_t radin, radout;
     Float_t me,alpha,yfact,ycell,bval;
     Float_t dE_IntBrehm, dE_ExtBrehmIn, dE_ExtBrehmOut, dE_Ionization;
     Float_t dE_Bsum;
     Float_t dE_IonizeIn, dE_IonizeOut; 
     std::vector<std::vector<Float_t> > Enumer;
     std::vector<TF1 *> fdistr;

     hamcEloss(const hamcEloss& eloss);
     hamcEloss& operator=(const hamcEloss& eloss);


#ifndef NODICT
ClassDef (hamcEloss, 0)   // electromagnetic radiation
#endif

};

#endif




   
