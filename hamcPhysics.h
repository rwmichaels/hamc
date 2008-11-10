#ifndef ROOT_hamcPhysics
#define ROOT_hamcPhysics

//  hamcPhysics   -- abstract base class for the physics of an experiment.
//  d. Jaunzeikare, R. Michaels  May 2008

#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;
// class hamcKine;   // see comments below
// class hamcRad;    //    "    "

class hamcPhysics {

  public:

     hamcPhysics();
     virtual ~hamcPhysics()=0;   // abstract class

     std::string GetName() const { return phy_name; };
     virtual Int_t Init(hamcExpt* expt);
     virtual Int_t Generate(hamcExpt* expt);    // Generate crsec and asy.
     virtual Float_t GetCrossSection() const { return crsec; };
     virtual Float_t GetAsymmetry() const { return asymmetry; };
     virtual Float_t GetdE_IntBrehm() const { return dE_IntBrehm; };
     virtual Float_t GetdE_ExtBrehm() const { return dE_ExtBrehm; };
     virtual Float_t GetdE_Ionization()  const { return dE_Ionization; };

  protected:

     virtual Int_t CrossSection();    
     virtual Int_t Asymmetry();
     virtual Int_t Radiate();

// You may want to create a separate kinematics class, as
// a member of hamcPhysics
//     hamcKine* kine;

// You may want to create a separate radiation class, as
// a member of hamcPhysics
//     hamcRad* radiation;

     Bool_t did_init;
     Float_t crsec, asymmetry;
     Float_t energy, theta_rad, theta_deg, phi_rad, phi_deg;
     Float_t dE_IntBrehm, dE_ExtBrehm, dE_Ionization;

     std::string phy_name;


  private: 


     hamcPhysics(const hamcPhysics& phys);
     hamcPhysics& operator=(const hamcPhysics& phys);


#ifndef NODICT
ClassDef (hamcPhysics, 0)   // physics
#endif

};

#endif



   
