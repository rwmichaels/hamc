#ifndef ROOT_hamcPhyPVDIS
#define ROOT_hamcPhyPVDIS

//  hamcPhyPVDIS   -- class for the PVDIS physics 

#include "Rtypes.h"
#include "hamcPhysics.h"
#include "TH1F.h"
#include <vector>
#include <string>
#include <map>



extern "C" { double cross_section__(double *, double *, double *, double *, double *);}
extern "C" {double r1998__(double *, double *);}
extern "C" {void getpdf_mrst2003c__(double *, double *, double *, double *);}
extern "C" {int NextUn__();}
extern "C" {void mrst2003c__(double *, double *, double *, double *,double *, double *, double *, double *,double *, double *, double *, double *);}
extern "C" {void readpdf_single__(char *, double *, double *, double *, double *,double *, double *, double *, double *,double *, double *, double *, double *, double *, double *, double *, double *,double *, double *);}


using namespace std;

class hamcExpt;

class hamcPhyPVDIS : public hamcPhysics {

  public:

     hamcPhyPVDIS();
     virtual ~hamcPhyPVDIS(); 

     Int_t Init(hamcExpt* expt);
     Int_t Generate(hamcExpt* expt);    // Generate crsec and asy.

  protected:

  private: 

     Bool_t didinit;
     double Z;        
     double A;


     hamcPhyPVDIS(const hamcPhyPVDIS& phys);
     hamcPhyPVDIS& operator=(const hamcPhyPVDIS& phys);


#ifndef NODICT
ClassDef (hamcPhyPVDIS, 0)   // PVDIS physics
#endif

};

#endif
