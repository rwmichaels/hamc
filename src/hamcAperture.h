#ifndef ROOT_hamcAperture
#define ROOT_hamcAperture

//  hamcAperture   -- aperture defining the acceptance
//                    for 1 location in spectrometer
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "TMath.h"
#include <vector>
#include <iostream>

class hamcAperture {  // an aperture that defines the acceptance
public:
    hamcAperture() : xcent(0),ycent(0) {};
    virtual ~hamcAperture() {};  
    virtual void Print() { 
      std::cout<<"Center " << xcent << "  "<<ycent<<std::endl;
    }
    void SetCenter(Float_t x, Float_t y) { xcent=x; ycent=y; };
    virtual Bool_t CheckAccept(Float_t x, Float_t y) const { return kTRUE; };
// An aperture may have material (degrader) or multiple materials.
    void DefineRadLen(Int_t idx, Float_t rl) { 
      if (idx < (Int_t)radlen.size()) {
         radlen[idx] = rl;
      } else {
	radlen.push_back(rl);
      }
    }
    virtual Float_t GetRadLen(Int_t idx) const { 
      if (idx < (Int_t)radlen.size()) return radlen[idx];
      return 0;
    }
    virtual Float_t GetRadLen(Float_t x, Float_t y) const {
 // default (virtual) method.  normally defined by inheritig class.
      if (radlen.size()>0) return radlen[0];
      return 0;
    }
protected:
    std::vector<Float_t> radlen;
    Float_t xcent,ycent;
#ifndef NODICT
ClassDef (hamcAperture, 0)   // Apertures used to define acceptance
#endif
};

class hamcBox : public hamcAperture {  // a box shaped aperture
public:
  hamcBox(Float_t xl, Float_t xh, Float_t yl, Float_t yh) : hamcAperture() {
    xlo = xl; xhi = xh, ylo = yl; yhi = yh;
  }
  ~hamcBox() { };
  void Print() { 
    std::cout<<"Aperture box "<<xlo<<"  "<<xhi<<"  "<<
       ylo<<"  "<<yhi<<std::endl;
    hamcAperture::Print();
  }
  Bool_t CheckAccept(Float_t x, Float_t y) const {
    Float_t xdif = x-xcent;
    Float_t ydif = y-ycent;
    return ( (xdif > xlo) && (xdif < xhi) && (ydif > ylo) && (ydif < yhi) );
  }
private:
  Float_t xlo, xhi, ylo, yhi;
#ifndef NODICT
ClassDef (hamcBox, 0)   // Box shaped aperture
#endif
};

class hamcCircle : public hamcAperture {  // a circular aperture
public:
  hamcCircle(Float_t rad) : hamcAperture() {
    radius_squared = rad*rad;
  }
  ~hamcCircle() { };
  void Print() { 
    std::cout<<"Aperture circle R = "<<TMath::Sqrt(radius_squared)<<std::endl;
    hamcAperture::Print();
  }
  Bool_t CheckAccept(Float_t x, Float_t y) const {
    Float_t xdif = x-xcent;
    Float_t ydif = y-ycent;
    return ( (xdif*xdif + ydif*ydif) < radius_squared );
  }
private:
  Float_t radius_squared;
#ifndef NODICT
ClassDef (hamcCircle, 0)   // Circle shaped aperture
#endif
};

class hamcTrapezoid : public hamcAperture { // a trapezoidal aperture (centered)
public:
  hamcTrapezoid(Float_t xlo, Float_t xhi,  Float_t ys, Float_t sl ) : hamcAperture() {
    xlow = xlo;  xhigh = xhi;   ysize = ys;  yslope = sl;
  }
  ~hamcTrapezoid() {};
  void Print() { 
    std::cout<<"Aperture trapezoid "<<xlow<<"  "<<xhigh;
    std::cout<<"  "<<ysize<<"  "<<yslope<<std::endl;
    hamcAperture::Print();
  }
  Bool_t CheckAccept(Float_t x, Float_t y) const {
  // No provision to subtract a center offset here.
    if (x < xlow || x > xhigh) return kFALSE;
    Float_t ychk = ysize + yslope*x;
    if (y < -1*ychk || y > ychk) return kFALSE;
    return kTRUE;
  }
private:
  Float_t xlow, xhigh, ysize, yslope;
#ifndef NODICT
ClassDef (hamcTrapezoid, 0)   // Trapezoid shaped aperture
#endif
};

class hamcPaulBox : public hamcAperture {  
// This is Paul's idea to have two boxes with different RL
public:
  hamcPaulBox(Float_t x1l, Float_t x1h, Float_t y1l, Float_t y1h, Float_t x2l, Float_t x2h, Float_t y2l, Float_t y2h) : hamcAperture() {
    xlo1 = x1l; xhi1 = x1h, ylo1 = y1l; yhi1 = y1h;
    xlo2 = x2l; xhi2 = x2h, ylo2 = y2l; yhi2 = y2h;
  }
  ~hamcPaulBox() { };
  void Print() { 
    std::cout<<"Paul's collimator with two boxes ";
    std::cout <<xlo1<<"  "<<xhi1<<"  "<<ylo1<<"  "<<yhi1<<std::endl;
    std::cout <<xlo2<<"  "<<xhi2<<"  "<<ylo2<<"  "<<yhi2<<std::endl;
    hamcAperture::Print();
  }
  Int_t WhichBox(Float_t x, Float_t y) const {
    Float_t xdif = x-xcent;
    Float_t ydif = y-ycent;
    Int_t idx = -1;
    if ((xdif > xlo1) && (xdif < xhi1) && 
        (ydif > ylo1) && (ydif < yhi1) ) idx=0;
    if ((xdif > xlo2) && (xdif < xhi2) && 
        (ydif > ylo2) && (ydif < yhi2) ) idx=1;
    return idx;
  }
  Float_t GetRadLen(Float_t x, Float_t y) const {
    return hamcAperture::GetRadLen(WhichBox(x,y));
  }
  Bool_t CheckAccept(Float_t x, Float_t y) const {
    Int_t idx = WhichBox(x,y);
    if (idx != -1) return kTRUE;
    return kFALSE;
  }
private:
  Float_t xlo1, xhi1, ylo1, yhi1;
  Float_t xlo2, xhi2, ylo2, yhi2;
#ifndef NODICT
ClassDef (hamcPaulBox, 0)   // Paul's complicated box.
#endif
};


#endif


   
