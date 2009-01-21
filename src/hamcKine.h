#ifndef ROOT_hamcKine
#define ROOT_hamcKine

//  hamcKine   -- Kinematics of the event
//  Events are generated uniform in phase space.
//  They may be weighted (later) by cross section.
//  R. Michaels  Nov 2008

#include "Rtypes.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <iostream>
#include <map>

class hamcExpt;

class hamcAccCell {  
// Utility class to account for how the acceptance is populated
public:
  hamcAccCell(Float_t thmin, Float_t thmax, Float_t phmin, Float_t phmax) {
    cmax = TMath::Cos(thmin);
    cmin = TMath::Cos(thmax);
    dcostheta = (cmax-cmin)/((Float_t)(numcell)); // >= 0
    pmin = phmin;
    pmax = phmax;
    dphi = (pmax-pmin)/((Float_t)(numcell));  // >= 0
    domega = dcostheta * dphi;
    xcnt = new Float_t[numcell*numcell];
    for (Int_t i=0; i<numcell*numcell; i++) xcnt[i] = 0;
  };
  void Increment(Float_t th, Float_t ph) {
    Int_t icell,jcell,ncell;
    Float_t cth = TMath::Cos(th);
    icell = ((Int_t)((cth - cmin)/dcostheta));
    jcell = ((Int_t)((ph-pmin)/dphi));
    ncell = icell*numcell + jcell;
    if (ncell >= 0 || ncell < numcell*numcell) xcnt[ncell] += 1.0;
  };
  Float_t GetCosTheta(Int_t icell) {
    return cmin + dcostheta*((Float_t)(icell));
  };
  Float_t GetTheta(Int_t icell) {
    Float_t cth = GetCosTheta(icell);
    return TMath::ACos(cth);
  };
  Float_t GetPhi(Int_t icell) {
    return pmin + dphi*((Float_t)(icell));
  };
  Float_t Num(Int_t icell) { 
    if (icell < 0 || icell >= numcell*numcell) return 0;
    return xcnt[icell];
  };
 ~hamcAccCell() { delete [] xcnt; };
 Float_t cmin, cmax, pmin, pmax, dcostheta, dphi, domega;
 static const Int_t numcell = 100;
 Float_t *xcnt;
};


class hamcKine {

  public:

     hamcKine();
     virtual ~hamcKine(); 

     Int_t Init(hamcExpt *expt);

     Int_t Init(std::string proc, 
          Float_t ebeam, Float_t theta_central, 
          Float_t mass_tgt,
	  Float_t thmin, Float_t thmax, 
          Float_t phmin, Float_t phmax,
          Float_t epmin=0, Float_t epmax=0);

     void SetDisDef(Float_t xlo, Float_t xhi, Float_t qsqlo, Float_t wsqlo);

     Int_t Generate(hamcExpt *expt);  // event generator
     Int_t Generate(Float_t ebeam, Float_t deafter); 
     Int_t IncrementAcceptance();
     void Print();

     hamcAccCell *acell;

// Event variables.
     Float_t energy, theta, phi;
     Float_t eprime, qsq, wsq, y, x, bigy;
     Float_t pprime, erecoil, dE_after;

  private:

     Int_t GenerateElastic();
     Int_t GenerateDis();
     Int_t ComputeKine(); 
     void Clear();
     Int_t CheckInit();

     Int_t iproc;
     static const Int_t proc_undef = -1 ;
     static const Int_t proc_elastic = 0;
     static const Int_t proc_dis = 1;
     static const Float_t mass_electron = 0.000511; // GeV
     static const Float_t mass_proton   = 0.938;    // GeV

     Bool_t did_init;
     Float_t ebeam, ebeam_central, theta_central;
     Float_t mass_tgt;  
     Float_t thmin, thmax, phmin, phmax, epmin, epmax;
     Float_t xbjlo, xbjhi, qsqlo, wsqlo;
     Float_t dP0_iter;

     hamcKine(const hamcKine& kine);
     hamcKine& operator=(const hamcKine& kine);


#ifndef NODICT
ClassDef (hamcKine, 0)   // kinematics class
#endif

};

#endif



   
