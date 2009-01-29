#ifndef ROOT_hamcAccAvg
#define ROOT_hamcAccAvg

//  hamcAccAvg   -- 

// R. Michaels  Jan 2009
// Class to account for how the acceptance is populated
// Averages over acceptance to obtain rate, asymmetry, solid angle.

#include "Rtypes.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <string>
#include <iostream>
#include <map>

class hamcExpt;

class hamcAccAvg {

  public:

     hamcAccAvg(Float_t thc, Float_t th1, Float_t th2, Float_t ph1, Float_t ph2);
     virtual ~hamcAccAvg(); 

     void Init();
     void InitHisto(Int_t i);
     void Increment(Float_t th, Float_t ph, Float_t rate, Float_t asy);
     void RunSummary();

     Float_t GetAsy();
     Float_t GetRate();
     Float_t GetOmega();
     Float_t GetTheta(Int_t icell);
     Float_t GetPhi(Int_t icell);
     Float_t GetNum(Int_t icell);


     void Print();

     static const Int_t numcell = 200; 
     static const Int_t ncellcut = 4; 

  private:

     Float_t central_angle, thmin, thmax, phmin, phmax, dtheta, dphi;
     Bool_t did_init, did_summary;
     Float_t *sumcnt, *sumrate, *sumasy;
     Float_t rate, asy, omega;
     Int_t debug;

     Int_t Check();

     TH2F *hacc;

     hamcAccAvg(const hamcAccAvg& acc);
     hamcAccAvg& operator=(const hamcAccAvg& acc);


#ifndef NODICT
ClassDef (hamcAccAvg, 0)   // Class to organize acceptance averaging
#endif

};

#endif



   
