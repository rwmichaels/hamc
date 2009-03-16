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
        for (Int_t ix = (Int_t)radlen.size(); ix < idx; ix++) {
  	    radlen.push_back(0);
	}
	radlen.push_back(rl);
      }
    }
    virtual Float_t GetRadLen(Int_t idx) const { 
      if (idx >= 0 && 
          idx < (Int_t)radlen.size()) return radlen[idx];
      return 0;
    }
    virtual Float_t GetRadLen(Float_t x, Float_t y) const {
 // default (virtual) method.  normally defined by inheriting class.
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
// The first one has a RL defined, the 2nd is vacuum normally.
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

class hamcPaulCollim : public hamcAperture {  
// This is Paul's idea to have an apertures a different RL
// The regular aperture is vacuum normally.
// The acceptance hugs the natural acceptance.
public:
  hamcPaulCollim(
    Float_t atx, Float_t aty, Float_t atrr, Float_t atcc, // A_T hole    
    Float_t r1, Float_t c1, // outer circle, main acceptance
    Float_t r2, Float_t c2, // inner circle, main acceptance
    Float_t ytp,            // y is below this line (on top) or above -1* this
    Float_t xrgt,           // largest pos. X (small angle side)
    Float_t x1, Float_t m1) // y is below this (top), above (bottom)
       : hamcAperture() {
    atxlo = atx; atyhi = aty;  atr = atrr;  atc = atcc;   
    rad1 = r1; ycent1 = c1;
    rad2 = r2; ycent2 = c2;
    vtop = ytp;  vright = xrgt;
    yline1 = x1; slope1 = m1;
  }
  ~hamcPaulCollim() { };
  void Print() { 
    std::cout<<"Paul's collimator with inner, outer radii"<<std::endl;
    std::cout<<"A_T hole :"<<std::endl;
    std::cout <<" (vert) X > "<<atxlo<<"  (hor) Y < "<<atyhi;
    std::cout <<"   rad "<<atr<<"  cent "<<atc<<std::endl;
    std::cout <<"outer circle "<<rad1<<"  "<<ycent1<<std::endl;
    std::cout <<"inner circle "<<rad2<<"  "<<ycent2<<std::endl;
    std::cout <<"Top/bottom +/= "<<vtop<<std::endl;
    std::cout <<"Right side cut "<<vright<<std::endl;    
    std::cout <<"Line (Champhor)  "<<yline1<<"  "<<slope1<<std::endl;
    std::cout <<"Rad lens "<<hamcAperture::GetRadLen(0);
    std::cout << "  "<<hamcAperture::GetRadLen(1);
    std::cout << "   size "<<radlen.size()<<std::endl;
    hamcAperture::Print();
  }
  Int_t WhichBox(Float_t x, Float_t y) const {
    Int_t debug=0;
    Float_t xdif = x-xcent;
    Float_t ydif = y-ycent;
    Float_t xtrial,ytrial,xtmp;
    if (debug) std::cout << "\n\n X, Y = "<<x<<"  "<<y<<std::endl;
    ytrial = ydif-atc;
    xtrial = atr*atr - ytrial*ytrial;
    if (xtrial > 0) {
      xtmp = TMath::Sqrt(xtrial);
      if ((xdif > atxlo) && (ydif < atyhi) && 
          (x < xtmp)) return 0;     // in A_T hole (a triangle with arc)
    } 
    // Now check if in main acceptance
    // Remember x is vertical (transport), y is horizontal
    ytrial = ydif-ycent1;
    xtrial = rad1*rad1 - ytrial*ytrial;
    if (debug) std::cout<<"chk1 "<<xtrial<<"  "<<x<<std::endl;
    if (xtrial < 0) return -1;   // outside outer circle
    xtrial = TMath::Sqrt(xtrial);
    if (debug) std::cout<<"chk2 "<<xtrial<<"  "<<x<<std::endl;
    if (x > xtrial || x < -1.0*xtrial) return -1;  // outside outer circle
    ytrial = ydif-ycent2;
    xtrial = rad2*rad2 - ytrial*ytrial;
    if (xtrial > 0) {  
      xtrial = TMath::Sqrt(xtrial);
      if (debug) std::cout<<"chk3 "<<xtrial<<"  "<<x<<std::endl;
      if ((x > 0 && x < xtrial) || 
          (x < 0 && x > -1.0*xtrial)) return -1; // outside innner circle
    }
    if (x > vtop) return -1;       // above upper line
    if (x < -1.0*vtop) return -1;  // beneath lower line
    if (y > vright) return -1;     // beyond right border
    xtrial = slope1*ydif + yline1; // Champhor line
    if (x > xtrial) return -1;     // above upper line
    if (x < -1.0*xtrial) return -1;// below lower line
    if (debug) std::cout<<"INSIDE "<<std::endl;
    return 1;  // inside main acceptance.
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
  Float_t x, atxlo, atyhi, atr, atc;
  Float_t rad1, ycent1, rad2, ycent2;
  Float_t yline1, slope1;
  Float_t vtop, vright;
#ifndef NODICT
ClassDef (hamcPaulBox, 0)   // Paul's semicircular, really complicated box.
#endif
};


#endif


   
