#ifndef ROOT_hamcSpecHRS
#define ROOT_hamcSpecHRS

//  hamcSpecHRS   -- class for the HRS spectometer
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcAperture.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>

// mnemonic constants
#define ITARGET       0
#define ITARGET_FULL 21
#define ICOLLIM       1
#define ICOLLIM2      2
#define ICOLLIM3     20
#define ISEPTIN       3
#define ISEPTOUT      4
#define IQ1EXIT       5
#define IDIPIN        6
#define IDIPEXIT      7
#define IQ3IN         8
#define IQ3EXIT       9
#define IFOCAL       10
#define IPLANE1      11
#define IPLANE2      12

#define LEFTHRS      1
#define RIGHTHRS     2


class hamcExpt;
class hamcTrack;
class hamcAperture;
class hamcTrans;
class hamcSpecHRS;

// class hamcDetector;

class hamcSpecBrk {
// Break points in spectrometer, where to perform transport
// and to possibly apply acceptance cuts and multiple scattering
public:
  hamcSpecBrk(Int_t wh, hamcAperture *app=0) {
    where = wh;
    aperture = app;
    accept_cut = kTRUE;
  }  
  Bool_t IsCut() { return accept_cut; }
  void DefineCut( Int_t icut ) { 
      std::cout << "hamcSpecBrk: DefineCut "<<where<<"  "<<icut<<std::endl;
      if (!icut) {
         accept_cut = kFALSE; 
      } else {
 	 accept_cut = kTRUE;
      }
  }
  void DefineMaterial(Int_t index, Float_t a, Float_t z, Float_t t) {
    if (aperture) aperture->DefineMaterial(index,a, z, t);
  }
  ~hamcSpecBrk() { 
    if (aperture) delete aperture;
  }
  void Print() {
    std::cout << "\n--------------------------------"<<std::endl;
    std::cout<<"Break point at "<<ItoName(where)<<std::endl;
    if (aperture) {
      aperture->Print();
    } else {
      std::cout << "no aperture"<<std::endl;
    }
    if ( !accept_cut ) std::cout<<"Acceptance cut turned off !"<<std::endl;
  }
  std::string ItoName(Int_t wh) {
      if (wh == ITARGET ) return "target";
      if (wh == ICOLLIM ) return "collimator";
      if (wh == ICOLLIM2) return "2nd collimator";
      if (wh == ICOLLIM3) return "empirical angle collimator";
      if (wh == ISEPTIN ) return "septum input";
      if (wh == ISEPTOUT) return "septum output";
      if (wh == IQ1EXIT ) return "Q2 exit ";
      if (wh == IDIPIN )  return "dipole input";
      if (wh == IDIPEXIT) return "dipole output";
      if (wh == IQ3IN  )  return "Q3 input";
      if (wh == IQ3EXIT)  return "Q3 output";
      if (wh == IFOCAL )  return "focal plane ";
      if (wh == IPLANE1)  return "plane 1";
      if (wh == IPLANE2)  return "plane 2";
      return "none";
  }
  Int_t where;
  Bool_t accept_cut;
  hamcAperture *aperture;
};

class hamcTransGuido;

class hamcSpecHRS {

  public:

     hamcSpecHRS(int id, Float_t pmom, Float_t angle);
     virtual ~hamcSpecHRS(); 

     Int_t Init(hamcExpt *expt);
     
     hamcAperture* Aperture(Int_t brk);
     hamcTrans* transport;
     hamcTransGuido *tguido;
     Bool_t IsMultScatt(Int_t brk);
     Float_t GetRadLen(Int_t brk);

     Int_t GetNumBrk();
     std::vector<hamcSpecBrk *> break_point;

     void Print();

// Initalization choices. 

     void UseCollimator();   // To use the front-end collmator 
     void UsePaulColl();     // To use Paul's 2-piece collimator.
     void UseAngleColl();     // To use empirically derived angle collimation

     void UseHRSOnly();      // HRS without septum
     void Use4degSeptum();   // To use with the 4 degree Septum
     void UseWarmSeptum();   // To use with the Warm Septum
     void UseColdSeptum();   // Use with Cold Septum (ca 2005)

     void UseMatrixTrans();  // To use 1st order TRANSPORT
     void UseLeroseTrans();  // To use LeRose Transfer functions
                             // for either HRS, warm or cold septum
     void UseGuidoTrans();   // To use Guido's parameterization.

// Query about setup

     Bool_t IsCollimated() { return collim_choice!=nocollim; };
     Bool_t IsPaulCollim() { return collim_choice==paul_coll; };
     Bool_t IsAngleCollim() { return collim_choice==angle_coll; };
     Bool_t IsColdSeptum()  { return sept_choice==coldsept; };
     Bool_t Is4degSeptum()  { return sept_choice==sept4deg; };
     Bool_t IsWarmSeptum()  { return sept_choice==warmsept; };
     Bool_t IsMatrixTrans() { return trans_choice==tmatrix; };
     Bool_t IsLeroseTrans() { return trans_choice==tlerose; };
     Bool_t IsGuidoTrans() { return use_guido; };

     Float_t GetP0() const { return P0;};
     Float_t GetP0Sigma() const { return P0_sigma; };
     Float_t GetScattAngle() const { return central_angle; };
     Float_t GetCollimDist() const { return collim_distance; };

     Int_t which_spectrom, numdet;
     std::string name, desc;
     Float_t P0, P0_sigma, central_angle;
     Float_t collim2_radlen1, collim2_radlen2;
     Float_t collim2_a, collim2_z, collim2_t;
 
  private:

// None yet:    std::vector<hamcDetector* > detectors;

     Float_t collim_distance;

     Int_t sept_choice, trans_choice, collim_choice;
     static const Int_t noseptum  =0;
     static const Int_t warmsept  =1;
     static const Int_t coldsept  =2;
     static const Int_t sept4deg  =3;
     static const Int_t tmatrix   =1;
     static const Int_t tlerose   =2;
     static const Int_t nocollim  =0;
     static const Int_t reg_coll  =1;
     static const Int_t paul_coll =2;
     static const Int_t angle_coll=3;
     Bool_t use_guido;

     Int_t BuildSpectrom();
     Int_t ChkIdx(Int_t idx);
     void AddBreakPoint(Int_t iwhere);

     hamcSpecHRS(const hamcSpecHRS& spec);
     hamcSpecHRS& operator=(const hamcSpecHRS& spec);

#ifndef NODICT
ClassDef (hamcSpecHRS, 0)   // HRS Spectrometer 
#endif


};


#endif


   
