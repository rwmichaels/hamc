#ifndef ROOT_hamcPhyWater
#define ROOT_hamcPhyWater

//  hamcPhyWater   -- class for the Water cell physics 

#include "Rtypes.h"
#include "hamcPhysics.h"
#include "TH1F.h"
#include <vector>
#include <string>
#include <map>

using namespace std;

class hamcExpt;

class hamcPhyWater : public hamcPhysics {

  public:

     hamcPhyWater();
     virtual ~hamcPhyWater(); 

     Int_t Init(hamcExpt* expt);
     Int_t Generate(hamcExpt* expt);    // Generate crsec and asy.

     Float_t O16CrossSection(Float_t energy, Float_t angle); 

     Float_t sig_elas_H(Float_t E_beam, Float_t theta, Float_t Q_sqr);
     Float_t GEp_dipol(Float_t Q_sqr);
     Float_t GMp_dipol(Float_t Q_sqr);
     Float_t GEn_dipol(Float_t Q_sqr);
     Float_t GMn_dipol(Float_t Q_sqr);
     Float_t GAp3(Float_t Q_sqr);
     Float_t GAp8(Float_t Q_sqr);
     Float_t GAp_dipol(Float_t Q_sqr);
     Float_t asym_H(Float_t theta, Float_t Q2);
     Int_t Drate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec);

  protected:

  private: 

     Bool_t didinit;
     static const Float_t MV2=0.71;
     static const Float_t lambda_n=5.6;
     static const Float_t M_A=1.001;
     static const Float_t gA=-1.2695;
     static const Float_t Ffact=0.463;
     static const Float_t Dfact=0.804;
     static const Float_t Alpha=0.007299;  // 1/137
     static const Float_t mypi=3.1415926;
     static const Float_t Mp=0.938;
     static const Float_t rhop=0.9881;
     static const Float_t lambda1u=-1.85E-5;
     static const Float_t lambda1d=3.705E-5;
     static const Float_t lambda2u=-0.0121;
     static const Float_t lambda2d=0.0026;
     static const Float_t hbc2=389.37966;  /* GeV2.mubarn */
     static const Float_t kappa=1.0300;
     static const Float_t kappap=1.0027;
     static const Float_t rho=1.0011;
     static const Float_t sw2=0.23117;  /* at Z0  mass in MSbar scheme */
     static const Float_t GFermi=1.16639E-05;         /* GeV-2 */
     static const Float_t mu_n=-1.9130428;
     static const Float_t mu_p=2.79284739;

     static const Int_t quick_check=1;

     hamcPhyWater(const hamcPhyWater& phys);
     hamcPhyWater& operator=(const hamcPhyWater& phys);


#ifndef NODICT
ClassDef (hamcPhyWater, 0)   // Water cell physics
#endif

};

#endif
