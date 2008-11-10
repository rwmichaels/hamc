#ifndef ROOT_hamcPhyPREX
#define ROOT_hamcPhyPREX

//  hamcPhyPREX   -- class for the PREX physics of an experiment.
//  D. Jaunzeikare, R. Michaels  May 2008
//  This is a work in progress !!  (Sept 2008)

#include "Rtypes.h"
#include "hamcPhysics.h"
#include "TH1F.h"
#include <vector>
#include <string>
#include <map>

using namespace std;

class hamcExpt;
//class hamcTarget;
//class hamcKine;
//class hamcRad;

class hamcPhyPREX : public hamcPhysics {

  public:

     hamcPhyPREX();
     virtual ~hamcPhyPREX(); 

     Int_t Init(hamcExpt* expt);
     Int_t Generate(hamcExpt* expt);    // Generate crsec and asy.

     /* Reads CrossSection and Asymmetry from Horowitz tables and saves in class variables*/
     Int_t CrossSection(Float_t energy, Float_t angle, Int_t stretch=0);    
     Int_t Asymmetry(Float_t energy, Float_t angle, Int_t stretch=0);

     Int_t Radiate();

     //moved to public because hamcExptPREX::Test needs to call it 

    /*Calculates CrossSection and Asymmetry using formulas*/
     Float_t CalculateCrossSection(Float_t energy, Float_t angle);
     Float_t CalculateAsymmetry(Float_t energy, Float_t angle);


     /*---- Radiative Corrections -----*/
     Int_t InitInternalRadCor();
     Int_t InitExternalRadCor();
     Float_t GenerateDeltaE_internal(); /*Get random deltaE from the elist */
     Float_t GenerateDeltaE_external();
     Int_t TestRadiate();

  protected:

  private: 
//     hamcExpt* gExpt; 

     Int_t LoadFiles();
     Int_t LoadHorowitzTable(vector<vector<Float_t> >&, vector<vector<Float_t> >&,Int_t stretch = 0);
     Int_t LoadFormFactorTable(); // qsq row, ffsq row
     vector<Float_t> qsq_row, ffsq_row;

     /*Helper methods*/
     Int_t FindAngleIndex(Float_t);
     Int_t FindEnergyIndex(Float_t);
     Float_t Interpolate(Float_t min, Float_t max, Float_t mid, Float_t val1, Float_t val2);
     Float_t CalculateQsq(Float_t energy, Float_t angle);

     Bool_t didinit;
//     Float_t testvar;   // a test variable

     hamcPhyPREX(const hamcPhyPREX& phys);
     hamcPhyPREX& operator=(const hamcPhyPREX& phys);

     /*3D vectors to save data from Horowitz tables*/
     vector< vector< vector<Float_t> > > crsc_tables;   
     vector< vector< vector<Float_t> > >asymmetry_tables;

     vector<Float_t> angle_row;  //to save the angle values of Horowitch table
     vector<Float_t> energy_row; //to save the energy values 


     /* We can probably move some of this to the base class, eventually */ 
    /*-----Radiative Corrections */
     Int_t InitRadCor(Float_t*, Int_t, Float_t, TH1F*); /*Divide the function (calculated by CalcualteIe() into cells and saves in elist */

     Float_t tlen;  // target length (RL)
     Float_t *elistInternal;
     Float_t *elistExternal;
         
     Int_t ncell1;
     Int_t ncell2;

     

#ifndef NODICT
ClassDef (hamcPhyPREX, 0)   // PREX physics
#endif

};

#endif
