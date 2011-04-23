#ifndef ROOT_hamcPhyPREX
#define ROOT_hamcPhyPREX

//  hamcPhyPREX   -- class for the PREX physics
//  D. Jaunzeikare, R. Michaels  May 2008

#include "Rtypes.h"
#include "hamcPhysics.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>

#define NOMODEL  0
#define HORPB    1
#define HORCA    2
#define HORSN    3
#define SI       4
#define NL3P06   5
#define SLY4     6
#define SIII     7
#define FSU      8
#define NL3      9
#define NL3M05  10

using namespace std;

class hamcExpt;

class hamcQuad {
// Quadrupole beamline element.
//         (x,x')    cm/rad  (not mrad)
//         (x',x)    rad/cm
 public:
  hamcQuad(Float_t bbp, Float_t rrad, Float_t llen) : 
// bbp = field (gauss) at pole,  rrad = radius of quad (cm), llen = length (cm)
    bp(bbp), rad(rrad), len(llen), focusx(1), pcut(0.01), huge(999999), debug(0) 
    {  grad = bp / rad;
      kappa0 = TMath::Sqrt(4.8e-10 * grad / 1.6e-3); }
  virtual ~hamcQuad() {};
  void SetFocusX() { focusx = 1; };
  void SetFocusY() { focusx = 0; };
  Float_t GetMatrix(Float_t pmom, Int_t index) {
    if (pmom < pcut) return huge;
    kappa = kappa0 / TMath::Sqrt(pmom);
    kL = kappa * len;
    exp1 = exp(kL);
    exp2 = exp(-1*kL);
    coshx = (exp1 + exp2)/2;
    sinhx = (exp1 - exp2)/2;
     if (debug) 
       std::cout << "focus "<<focusx<<"  pmom "<<pmom<<"  kappa "<<kappa<<"  kL "<<kL<<std::endl;
    if (focusx) {
      if (index == 0) return TMath::Cos(kL);
      if (index == 1) return (1/kappa)*TMath::Sin(kL);
      if (index == 2) return -1*kappa*TMath::Sin(kL);
      if (index == 3) return TMath::Cos(kL);
      if (index == 4) return coshx;
      if (index == 5) return (1/kappa)*sinhx;
      if (index == 6) return kappa*sinhx;
      if (index == 7) return coshx;
    } else {
      if (index == 0) return coshx;
      if (index == 1) return (1/kappa)*sinhx;
      if (index == 2) return kappa*sinhx;
      if (index == 3) return coshx;
      if (index == 4) return TMath::Cos(kL);
      if (index == 5) return (1/kappa)*TMath::Sin(kL);
      if (index == 6) return -1*kappa*TMath::Sin(kL);
      if (index == 7) return TMath::Cos(kL);
    }
    return 0;
  }
  void Print() {
    std::cout << "Quad values "<<std::endl;
    std::cout << "grad "<<grad<<"    kappa "<<kappa<<std::endl;
    if (focusx) {
      std::cout<<"focus in X"<<std::endl;
    } else {
      std::cout<<"focus in Y"<<std::endl;
    }
  }
 private:
  Float_t bp, rad, len;
  Int_t focusx;
  Float_t pcut, huge;
  Int_t debug;
  Float_t grad, kappa0;
  Float_t exp1, exp2;
  Float_t kappa, kL, coshx, sinhx;
#ifndef NODICT
ClassDef (hamcQuad, 0)   // Quadrupole beamline element
#endif
};

class hamcPhyPREX : public hamcPhysics {

  public:

     hamcPhyPREX();
     virtual ~hamcPhyPREX(); 

     Int_t Init(hamcExpt* expt);
     Int_t Generate(hamcExpt* expt);    // Generate crsec and asy.

// The index "i" may point to a model or parameter set.
     Float_t GetCrossSection(Int_t i=0) const;
     Float_t GetAsymmetry(Int_t i=0) const;
     /* Reads CrossSection and Asymmetry from Horowitz tables and saves in class variables*/
     Int_t CrossSection(Float_t energy, Float_t angle, Int_t stretch=0);    
     Int_t Asymmetry(Float_t energy, Float_t angle, Int_t stretch=0);
     Int_t Drate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec);

    /*Calculates CrossSection and Asymmetry using formulas*/
     Float_t CalculateCrossSection(Int_t nuc, Float_t energy, Float_t angle);
     Float_t CalculateAsymmetry(Int_t nuc);  // Born asymmetry.
     Float_t CalculateDrate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec);
     Int_t SetModel(Int_t modeln);  // To set the model used

  protected:

  private: 

     Int_t LoadFiles();
     Int_t LoadHorowitzTable(vector<vector<Float_t> >&, vector<vector<Float_t> >&,Int_t stretch = 0);
     Int_t LoadFormFactorTable(); // qsq row, ffsq row
     vector<Float_t> qsq_row, ffsq_row;

     Int_t LoadC12FormFactorTable(); // qsq row, ffsq row
     vector<Float_t> c12qsq_row, c12ffsq_row;  // same for C12

     /*Helper methods*/
     Int_t FindAngleIndex(Float_t);
     Int_t FindEnergyIndex(Float_t);
     Float_t Interpolate(Float_t min, Float_t max, Float_t mid, Float_t val1, Float_t val2);
     Float_t InterpAsym(Float_t ene, Float_t angle_rad);
     void PrintAsymTable();

     hamcQuad *quad1, *quad2;

     Float_t qsq;
     Float_t asy0, asy1;  // unstretched and streteched R_N asymmetries 
     
     Float_t drate;
     Bool_t didWriteAngle;

     // Now considering 8 models.
     // The first one is "old" (~10 years)

     // Model        Rn(fm)    A(850 MeV, 6 deg)  A(1050 MeV, 5deg)
     // 1. horpb.dat and horpb1.dat
     // 2. SI          5.492   .731               .753
     // 3. NL3p06      5.604   .700               .721
     // 4. SLY4        5.618   .698               .719
     // 5. SIII        5.646   .683               .702
     // 6. FSU         5.679   .677               .695
     // 7. NL3         5.740   .655               .672
     // 8. NL3m05      5.851   .617               .631
     
     Int_t whichmodel;

     Bool_t didinit;

     hamcPhyPREX(const hamcPhyPREX& phys);
     hamcPhyPREX& operator=(const hamcPhyPREX& phys);

     /*3D vectors to save data from Horowitz tables*/
     vector< vector< vector<Float_t> > > crsc_tables;   
     vector< vector< vector<Float_t> > >asymmetry_tables;

     vector<Float_t> angle_row;  //to save the angle values of Horowitch table
     vector<Float_t> energy_row; //to save the energy values 

     TH1F *histphi; 
     TH1F *histene;
     TH1F *histprob;
     TH2F *histrast;
     TH1F *histinacc;
     TH1F *histz1a, *histz1b, *histz1c, *histz1d, *histz1e, *histz1f;
     TH1F *histz2a, *histz2b, *histz2c, *histz2d, *histz2e, *histz2f;
     TH1F *histz3a, *histz3b, *histz3c, *histz3d, *histz3e, *histz3f;
     TH1F *histz4a, *histz4b, *histz4c, *histz4d, *histz4e, *histz4f;
     TH1F *histz5a, *histz5b;
     TH1F *hpph1,*hpph2,*hpph3,*hpph4,*hpph6,*hd1;
     TH1F *hfom1,*hfom2,*hfom3,*hfom4,*hfom5, *hfom6;
     TH1F *hfom7, *hfom8, *hfom9, *hfom10, *hfom11, *hfom12, *hfom13;
     TH2F *hpph5;
     static const Int_t histo_test=1; // to test(1) or not(0) this code
     static const Int_t accept_test=0; // to test(1) or not(0) this acceptance
     static const Int_t quick_feasibility=0; 
     static const Int_t quick_check=0; 
     static const Int_t quick_fom=0; 
     static const Int_t power_integ=0; 
     static const Int_t neutron_power=0; 

#ifndef NODICT
ClassDef (hamcPhyPREX, 0)   // PREX physics
#endif

};

#endif
