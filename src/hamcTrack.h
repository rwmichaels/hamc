#ifndef ROOT_hamcTrack
#define ROOT_hamcTrack

//  hamcTrack   -- base class for a track
//  Track is defined by TRANSPORT coordinates and by its 4-vector
//  in lab frame.  Note, the 4-vector is only valid in the field-free
//  region where it originates -- at the target.  In contrast,
//  the TRANSPORT vector may be modified by a hamcTrans object.
//
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcExpt.h"
#include <string>
#include <iostream>
#include <vector>

class hamcExpt;
class hamcTrans;
class hamcAperture;
class hamcBeam;
class hamcSpectrom;

class hamcTransVect {
// Transport vector.  Utility class needed by hamcTrack
public:
  hamcTransVect() { Load(0.,0.,0.,0.,0.); };
  hamcTransVect(Float_t x, Float_t th, Float_t y, Float_t ph, Float_t dp) {
      Load(x,th,y,ph,dp);
  }
  std::vector<Float_t> Get() { return vtrans; };
  Float_t GetX() const { return vtrans[0]; };
  Float_t GetTheta() const { return vtrans[1]; };
  Float_t GetY() const { return vtrans[2]; };
  Float_t GetPhi() const { return vtrans[3]; };
  Float_t GetDpp() const { return vtrans[4]; };
  Float_t GetZ() const { return zscatt; };
  void PutX(Float_t x) { vtrans[0]=x; };
  void PutTheta(Float_t th) { vtrans[1]=th; };
  void PutY(Float_t y) { vtrans[2]=y; };
  void PutZ(Float_t z) { zscatt = z; };
  void PutPhi(Float_t ph) { vtrans[3]=ph; };
  void PutDpp(Float_t dp) { vtrans[4]=dp; };
  void Load(Float_t *v) { Load(v[0],v[1],v[2],v[3],v[4]); };
  void Load(Float_t x, Float_t th, Float_t y, Float_t ph, Float_t dp) {
    vtrans.clear();
    vtrans.push_back(x);
    vtrans.push_back(th);
    vtrans.push_back(y);
    vtrans.push_back(ph);
    vtrans.push_back(dp);
  };
  void Clear() { Load(0.,0.,0.,0.,0.); };
  void Load(std::vector<Float_t> data) {
    vtrans.clear();
    for (Int_t i=0; i<5; i++) vtrans.push_back(data[i]);
  };
  void AddToTheta(Float_t th) { vtrans[1] += th; };
  void AddToPhi(Float_t ph) { vtrans[3] += ph; };
  void Print() {
    std::cout << "Transport vector. ";
    std::cout << "  Size :  "<<vtrans.size()<<std::endl;
    std::cout << "X = "<<GetX()<<"   Theta = "<<GetTheta();
    std::cout << "   Y = "<<GetY()<<"   Phi = "<<GetPhi();
    std::cout << "   dpp = "<<GetDpp()<<std::endl;
    std::cout << "Z location (on beamline) "<<zscatt<<std::endl;
  };
  hamcTransVect& operator=(const hamcTransVect& rhs) {
    vtrans.clear();
    vtrans.push_back(rhs.GetX());
    vtrans.push_back(rhs.GetTheta());
    vtrans.push_back(rhs.GetY());
    vtrans.push_back(rhs.GetPhi());
    vtrans.push_back(rhs.GetDpp());
    return *this;
  };  
  Float_t zscatt;  // Z loc of scatt; only meaningful at target
private:
  std::vector<Float_t> vtrans;
  hamcTransVect(const hamcTransVect& tv);
};


class hamcTrack {
// Base class for a track

  public:

     hamcTrack(std::string pid);
     hamcTrack(std::string pid, Float_t energy, Float_t x, Float_t theta, Float_t y, Float_t phi, Float_t dp);
     virtual ~hamcTrack()=0;   // abstract

     virtual void Print();

// Transport vectors:  present and origin.
     hamcTransVect *tvect, *tvect_orig;
// The following is a copy of transport vector "tvect" for present event
     Float_t xtrans, thtrans, ytrans, phtrans, dpptrans, ztrans;
// The following transport right after scattering.
     Float_t x0,th0,y0,ph0,dpp0,z0;

// Index for track origin
     Int_t origin;     // ITARGET, ICOLLIM, etc.(def'n in hamcSpecHRS.h)
                           
// Energy and momentum in lab frame (at the target only)
     virtual Float_t GetEnergy() const { return energy; };
     virtual Float_t GetE0() const { return E0; };
     virtual Float_t GetPmom() const { return pmom; };
     virtual Float_t GetPx() const { return plab_x; };
     virtual Float_t GetPy() const { return plab_y; };
     virtual Float_t GetPz() const { return plab_z; };

// Transport coordinates
     virtual Float_t GetTransX() const { return tvect->GetX(); };
     virtual Float_t GetTransTheta() const { return tvect->GetTheta(); };
     virtual Float_t GetTransY() const { return tvect->GetY(); };
     virtual Float_t GetTransPhi() const { return tvect->GetPhi(); };
     virtual Float_t GetTransDp() const { return tvect->GetDpp(); };
     virtual Float_t GetTransZ() const { return tvect->GetZ(); };

// Polar coordinates in Lab (Z = beam axis)
     virtual Float_t GetTheta() const { return theta; }; // rad
     virtual Float_t GetPhi() const { return phi; }; //rad
     virtual Float_t GetThetaDeg() const { return 180*theta/PI; }; // degr
     virtual Float_t GetPhiDeg() const { return 180*phi/PI; }; //degr

// TRANSPORT.  This modifies tvect only.
     virtual Int_t Transport(const hamcTrans *trans, Int_t where);

// To check if track is in the acceptance at 'where'
     virtual Bool_t InAccept(const hamcAperture *acc);

// Multiple scattering
     virtual void MultScatt(const hamcExpt *expt, Int_t where);
     virtual void MultScatt(const hamcAperture *app, Int_t where);
     virtual void MultScatt(Float_t radlen, Int_t where);


  protected:

     virtual Int_t Init();
     virtual Int_t InitMass();
     virtual void UpdateTrans();

     std::string trktype;  // "beam", "out", etc
     std::string pid;      // "electron", "proton", etc

     Float_t E0, P0;       // central energy, momentum
     Float_t E0sigma;      // spread in energy
     Float_t P0sigma;      // spread in momentum
     Float_t energy, pmom; // energy, momentum 
     Float_t mass;         // mass 
     Float_t theta, phi;   // polar angles
     Float_t dP0_iter;     // offset in P0 (% of P0, for iteration)
     Float_t plab_x, plab_y, plab_z;   // Pmom components in lab
     Float_t xdet, ydet;   // detector frame vars. (in focal plane)
     Int_t ms_collim;      // mult. scatt. in collim ? (1/0)

     Bool_t did_init, inaccept;     
      
     static const Int_t debug=0;


  private: 

     hamcTrack& operator=(const hamcTrack& trk);
     hamcTrack(const hamcTrack& trk);

#ifndef NODICT
ClassDef (hamcTrack, 0)   // Track
#endif

};

#endif


   
