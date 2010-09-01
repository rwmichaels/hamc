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
#include "TH1F.h"
#include <map>
#include "hamcBeam.h"
#include "hamcTrackOut.h"

class hamcExpt;

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
     Int_t Generate(hamcExpt *expt);     // event generator
     Int_t GenerateOut(hamcExpt *expt);  // modify track for MS and Eloss, update Qsq.
     Int_t Generate(Float_t ebeam, Float_t deafter); 
     void Print();

// Event variables.
     Float_t energy, theta, phi;
     Float_t eprime, qsq, wsq, y, x, bigy;
     Float_t qsq_obs, qsq_atrk, qsqfr;
     Float_t pprime, erecoil, dE_int;
     Float_t theta_trans, phi_trans;
     Float_t scat_ang;
     hamcBeam *beam;
     hamcTrackOut *track;
     TH1F *eprime_gen;

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

     Float_t iteration;

     hamcKine(const hamcKine& kine);
     hamcKine& operator=(const hamcKine& kine);


#ifndef NODICT
ClassDef (hamcKine, 0)   // kinematics class
#endif

};

#endif



   
