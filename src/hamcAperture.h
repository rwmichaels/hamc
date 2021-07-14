#ifndef ROOT_hamcAperture
#define ROOT_hamcAperture

//  hamcAperture   -- aperture defining the acceptance
//                    for 1 location in spectrometer
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include "hamcTrack.h"
#define LEFTHRS 1

using namespace std;

class hamcAperture {  // an aperture that defines the acceptance
public:
    hamcAperture() : xcent(0),ycent(0) {};
    virtual ~hamcAperture() {};  
    virtual void Print() const { 
      std::cout<<"Center " << xcent << "  "<<ycent<<std::endl;
    }
    void SetCenter(Float_t x, Float_t y) { xcent=x; ycent=y; };
    virtual Bool_t CheckAccept(Float_t x, Float_t y) const { return kTRUE; };
    virtual Bool_t CheckAccept(hamcTrack *trk) const { 
      if (!trk) return kFALSE;
      return CheckAccept(trk->tvect->GetX(), trk->tvect->GetY());
    };
// An aperture may have material (degrader) or multiple materials.
    void DefineMaterial(Int_t idx, Float_t ta, Float_t tz, Float_t tt) { 
      if (idx < (Int_t) a.size()) {
         a[idx] = ta;
         z[idx] = tz;
         t[idx] = tt;
      } else {
        for (Int_t ix = (Int_t)a.size(); ix < idx; ix++) {
         a.push_back(0);
         z.push_back(0);
         t.push_back(0);
	}
	a.push_back(ta);
	z.push_back(tz);
	t.push_back(tt);
      }

    }
    virtual Float_t GetRadLen(Int_t idx) const { 
	double X0;
      if (idx >= 0 && 
          idx < (Int_t)a.size()){
	  X0 = 716.4*a[idx]/(z[idx]*(z[idx]+1.0)*log(287.0/sqrt(z[idx])));
	  return t[idx]/X0;
      }
      return 0;
    }
    virtual Float_t GetA(Int_t idx) const { 
      if (idx >= 0 && 
          idx < (Int_t)a.size()){
	  return a[idx];
      }
      return 0;
    }
    virtual Float_t GetZ(Int_t idx) const { 
      if (idx >= 0 && 
          idx < (Int_t)z.size()){
	  return z[idx];
      }
      return 0;
    }
    virtual Float_t Gett(Int_t idx) const { 
      if (idx >= 0 && 
          idx < (Int_t)t.size()){
	  return t[idx];
      }
      return 0;
    }

    virtual Float_t GetRadLen(Float_t x, Float_t y) const {
	double X0;
      if (a.size()>0){
        X0 = 716.4*a[0]/(z[0]*(z[0]+1.0)*log(287.0/sqrt(z[0])));
	  return t[0]/X0;
      }
      return 0;
    }

    virtual Float_t GetA(Float_t x, Float_t y) const {
 // default (virtual) method.  normally defined by inheriting class.
      if (a.size()>0){
        return a[0];
      }
      return 0;
    }
    virtual Float_t GetZ(Float_t x, Float_t y) const {
 // default (virtual) method.  normally defined by inheriting class.
      if (z.size()>0){
        return z[0];
      }
      return 0;
    }
    virtual Float_t Gett(Float_t x, Float_t y) const {
 // default (virtual) method.  normally defined by inheriting class.
      if (t.size()>0){
        return t[0];
      }
      return 0;
    }

    virtual void FlipPhi() { };
protected:
    std::vector<Float_t> a;
    std::vector<Float_t> z;
    std::vector<Float_t> t;
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
  void Print() const { 
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
  void Print() const { 
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
  void Print() const { 
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
  void Print() const { 
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
  Float_t GetA(Float_t x, Float_t y) const {
    return hamcAperture::GetA(WhichBox(x,y));
  }
  Float_t GetZ(Float_t x, Float_t y) const {
    return hamcAperture::GetZ(WhichBox(x,y));
  }
  Float_t Gett(Float_t x, Float_t y) const {
    return hamcAperture::Gett(WhichBox(x,y));
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
    ysign = 1.0; 
 }
  ~hamcPaulCollim() { };
  void Print() const { 
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
//    std::cout << "   size "<<radlen.size()<<std::endl;
    std::cout << "   size "<<a.size()<<std::endl;
    hamcAperture::Print();
  }
  Int_t WhichBox(Float_t x, Float_t y) const {
    Int_t debug=0;
    Float_t xdif = x-xcent;
    Float_t ydif = ysign*(y-ycent);  
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
    if (ydif > vright) return -1;     // beyond right border
    xtrial = slope1*ydif + yline1; // Champhor line
    if (x > xtrial) return -1;     // above upper line
    if (x < -1.0*xtrial) return -1;// below lower line
    if (debug) std::cout<<"INSIDE "<<std::endl;
    return 1;  // inside main acceptance.
  }
  Float_t GetRadLen(Float_t x, Float_t y) const {
    return hamcAperture::GetRadLen(WhichBox(x,y));
  }
  Float_t GetA(Float_t x, Float_t y) const {
    return hamcAperture::GetA(WhichBox(x,y));
  }
  Float_t GetZ(Float_t x, Float_t y) const {
    return hamcAperture::GetZ(WhichBox(x,y));
  }
  Float_t Gett(Float_t x, Float_t y) const {
    return hamcAperture::Gett(WhichBox(x,y));
  }

  Bool_t CheckAccept(Float_t x, Float_t y) const {
    Int_t idx = WhichBox(x,y);
    if (idx != -1) return kTRUE;
    return kFALSE;
  }
  virtual void FlipPhi() { 
      ysign = -1;  
  };
private:
  Float_t x, atxlo, atyhi, atr, atc;
  Float_t rad1, ycent1, rad2, ycent2;
  Float_t yline1, slope1;
  Float_t vtop, vright;
  Float_t ysign;
#ifndef NODICT
ClassDef (hamcPaulBox, 0)   // Paul's semicircular, really complicated box.
#endif
  };

  class hamcPREX2Coll : public hamcAperture {  
    // The PREX-II collimator
  public:
    hamcPREX2Coll(Int_t which_spectrom) : hamcAperture() {  
        rmin = 0.175;
        rmax = 0.315;
        rmin_c = 0.147 + 0.352*rmin;
        rmax_c = 0.147 + 0.352*rmax;
        chamfer1_m = 1.88;
        chamfer1_b = 0.0722;
        chamfer1_s = 0.04;
        chamfer2_m = 10.;
        chamfer2_b = 0.05;
        chamfer2_s = 0.03456;
        xmax = 0.089;
        ymin = -0.03975;
        osign=-1;
        if(which_spectrom == LEFTHRS) {
// actually the same transport model is used for R-HRS and L-HRS,
// so no flipping here; only the sign of tg_ph is flipped for comparison 
// to real data at the target.
	  osign = -1.0; 
        }
    }
    ~hamcPREX2Coll() {};
    Int_t WhichBox( Float_t x, Float_t y ) const {
       if( fabs(x) > xmax) return -1;
       if( osign*y < ymin ) return -1;
/*
       if( x < ((-chamfer1_m*(y+chamfer1_s) - chamfer1_b)) ) return -1;
       if( x > ((chamfer1_m*(y+chamfer1_s) + chamfer1_b)) ) return -1;

       if( x < ((-chamfer2_m*(y+chamfer2_s) - chamfer2_b)) ) return -1;
       if( x > ((chamfer2_m*(y+chamfer2_s) + chamfer2_b)) ) return -1;

       if( x*x + (y+rmin_c)*(y+rmin_c) < rmin*rmin ) return -1;
       if( x*x + (y+rmax_c)*(y+rmax_c) > rmax*rmax ) return -1;
*/
       if( x < ((-chamfer1_m*(osign*y+chamfer1_s) - chamfer1_b)) ) return -1;//New
       if( x > ((chamfer1_m*(osign*y+chamfer1_s) + chamfer1_b)) ) return -1;//New

       if( x < ((-chamfer2_m*(osign*y+chamfer2_s) - chamfer2_b)) ) return -1;//New
       if( x > ((chamfer2_m*(osign*y+chamfer2_s) + chamfer2_b)) ) return -1;//New

       if( sqrt(x*x + (osign*y+rmin_c)*(osign*y+rmin_c)) < rmin ) return -1;//New
       if( sqrt(x*x + (osign*y+rmax_c)*(osign*y+rmax_c)) > rmax ) return -1;//New

       return 1;
    };
    Bool_t CheckAccept(Float_t x, Float_t y) const {
      Int_t idx = WhichBox(x,y);
      if (idx != -1) return kTRUE;
      return kFALSE;
    }
    private:
      Float_t osign; // -1 for L-HRS, +1 for R-HRS
      Float_t rmin,rmax,rmin_c,rmax_c;
      Float_t chamfer1_m,chamfer1_b,chamfer1_s;
      Float_t chamfer2_m,chamfer2_b,chamfer2_s;
      Float_t xmax,ymin;
#ifndef NODICT
ClassDef (hamcPREX2Coll, 0)   // PREXII collimator
#endif
  };

class hamcAngleCollim1 : public hamcAperture {  

// this version is for PREX-1

// This is a purely empirical collimation for the HRS.
// The target angles tg_th and tg_ph (called th0 and phi0 in hamc)
// are checked to see if they are inside or outside a polygon
// defined by points hpt[],vpt[] where "h" means the horizontal
// angle (tg_ph) and "v" means vertical.
public:
  hamcAngleCollim1() : hamcAperture() {
    npts = 0;
    blowup = 1.0;   
    infty = 999999;
    debug = 0;
    Init();
  }
  ~hamcAngleCollim1() { };
  void Init() {
    std::cout <<" Into AngleCollim1 init"<<std::endl;
    // These points must be in the order of the polygon (going CCW for example).
    // Units are radians.  "v" is vertical.

#ifdef TEST1
    hpt.push_back(0);   vpt.push_back(-0.01);  
    hpt.push_back(0);   vpt.push_back(0.01);   
    hpt.push_back(0.005);   vpt.push_back(0.03);
    hpt.push_back(0.02);   vpt.push_back(0.01); 
    hpt.push_back(0.02);   vpt.push_back(-0.01);
    hpt.push_back(0.005);   vpt.push_back(-0.03);
#endif

    // This is the "golden" version for thin C12
    // Modified (Sept 17) to allow smaller angles

    hpt.push_back(-0.013);  vpt.push_back(0.);    
    hpt.push_back(-0.013);   vpt.push_back(0.007); 
    hpt.push_back(-0.013);  vpt.push_back(0.0165);
    hpt.push_back(-0.013); vpt.push_back(0.023);  
    hpt.push_back(-0.013);  vpt.push_back(0.029);  
    hpt.push_back(-0.0114); vpt.push_back(0.035);  
    hpt.push_back(-0.0081); vpt.push_back(0.0405); 
    hpt.push_back(-0.004);  vpt.push_back(0.048);  
    hpt.push_back(0);       vpt.push_back(0.049);  
    hpt.push_back(0.005);   vpt.push_back(0.048);  
    hpt.push_back(0.01);    vpt.push_back(0.049);  
    hpt.push_back(0.0141);  vpt.push_back(0.0405); 
    hpt.push_back(0.0182);  vpt.push_back(0.031);  
    hpt.push_back(0.0208);  vpt.push_back(0.0275); 
    hpt.push_back(0.0240);  vpt.push_back(0.017);  
    hpt.push_back(0.0258);  vpt.push_back(0.008);  
    hpt.push_back(0.0266);  vpt.push_back(0.0);    
    hpt.push_back(0.0256);  vpt.push_back(-0.0095);
    hpt.push_back(0.024);   vpt.push_back(-0.0196);
    hpt.push_back(0.0209);  vpt.push_back(-0.030); 
    hpt.push_back(0.0182);  vpt.push_back(-0.040); 
    hpt.push_back(0.0142);  vpt.push_back(-0.0495);
    hpt.push_back(0.010);   vpt.push_back(-0.050); 
    hpt.push_back(0.005);   vpt.push_back(-0.049); 
    hpt.push_back(0.0);     vpt.push_back(-0.050); 
    hpt.push_back(-0.004);  vpt.push_back(-0.049); 
    hpt.push_back(-0.008);  vpt.push_back(-0.040); 
    hpt.push_back(-0.0113);  vpt.push_back(-0.037);
    hpt.push_back(-0.013);  vpt.push_back(-0.032);
    hpt.push_back(-0.013);  vpt.push_back(-0.023);
    hpt.push_back(-0.013);   vpt.push_back(-0.013);

   
    npts = vpt.size();
  
    Float_t osign = -1;  // +1 for L-HRS, -1 for R-HRS
    ohmin = 99999;
    ohmax = -99999;
    ovmin = 99999;
    ovmax = -99999;

    Float_t hshift = 0.0;

    for (Int_t i = 0; i < npts; i++) {
  // Possible Sign flip for an HRS (compared to hamc)
        hpt[i] = osign*hpt[i] + hshift;
        hpt[i] = blowup * hpt[i];
        vpt[i] = blowup * vpt[i];
        if (hpt[i] < ohmin) ohmin = hpt[i];
        if (hpt[i] > ohmax) ohmax = hpt[i];
        if (vpt[i] < ovmin) ovmin = vpt[i];
        if (vpt[i] > ovmax) ovmax = vpt[i];

	//	std::cout << "Emp. angle "<<i<<"  "<<hpt[i]<<"  "<<vpt[i]<<"  signs "<<hsign[i]<<"   "<<vsign[i]<<std::endl;
    }

     Print();
  }
  void Print() const { 
    std::cout<<"Empirical angle collimation. **** Num points = "<<npts<<std::endl;
    for (Int_t i = 0; i < npts; i++) {
      std::cout << "Pt "<<i<<"  "<<hpt[i]<<"  "<<vpt[i];
     }
    hamcAperture::Print();
  }
  Int_t WhichBox(Float_t vert, Float_t horiz) const {
 // Check angles from target.  Units are radians.
 // Here, horiz is horizontal angle, vert = vertical

    if (debug) std::cout << "Into WhichBox vert,horiz =  "<<vert<<"  "<<horiz<<std::endl;
    if (horiz < ohmin || horiz > ohmax) return -1;
    if (vert < ovmin || vert > ovmax) return -1;

    Int_t i, j, c = 0;
    for (i = 0, j = npts-1; i < npts; j = i++) {
      if ( ((vpt[i]>vert) != (vpt[j]>vert)) &&
	   (horiz < (hpt[j]-hpt[i]) * (vert-vpt[i]) / (vpt[j]-vpt[i]) + hpt[i]) )
	c = !c;
    }
    if (c == 0) return -1;
    return 1;

  }
  Bool_t CheckAccept(hamcTrack *trk) const {
    if (!trk) return kFALSE;
    Float_t theta = trk->tvect->GetTheta();
    Float_t phi = trk->tvect->GetPhi();
    Int_t idx = WhichBox(theta,phi);
    if (idx == 1) return kTRUE;
    return kFALSE;
  }
private:
  Int_t npts;
  Int_t debug;
  Float_t blowup, infty;
  Float_t ohmin, ohmax; 
  Float_t ovmin, ovmax; 
  std::vector<Float_t> vpt,hpt;

#ifndef NODICT
ClassDef (hamcAngleCollim1, 0)   // Empirical angle collimation.
#endif
};

class hamcAngleCollim : public hamcAperture {  
  //
  // Updated for prex-2
  //
// This is a purely empirical collimation for the HRS.
// The target angles tg_th and tg_ph (called th0 and phi0 in hamc)
// are checked to see if they are inside or outside a polygon
// defined by points hpt[],vpt[] where "h" means the horizontal
// angle (tg_ph) and "v" means vertical.
public:
  hamcAngleCollim() : hamcAperture() {
    npts = 0;
    blowup = 1.0;   
    infty = 999999;
    debug = 0;
    Init();
  }
  ~hamcAngleCollim() { };
  void Init() {
    std::cout <<" Into AngleCollim init"<<std::endl;
    // These points must be in the order of the polygon (going CCW for example).
    // Units are radians.  "v" is vertical.

#ifdef TEST1
    hpt.push_back(0);   vpt.push_back(-0.01);  
    hpt.push_back(0);   vpt.push_back(0.01);   
    hpt.push_back(0.005);   vpt.push_back(0.03);
    hpt.push_back(0.02);   vpt.push_back(0.01); 
    hpt.push_back(0.02);   vpt.push_back(-0.01);
    hpt.push_back(0.005);   vpt.push_back(-0.03);
#endif

    // A first try at PREX-2 collimation.  Sept 24, 2020
    //based on PREX-2 old database
/*    hpt.push_back(0.01300);  vpt.push_back(0.04250);
    hpt.push_back(0.01200);  vpt.push_back(0.04250);
    hpt.push_back(0.01100);  vpt.push_back(0.04250);
    hpt.push_back(0.01000);  vpt.push_back(0.04250);
    hpt.push_back(0.00900);  vpt.push_back(0.04250);
    hpt.push_back(0.00800);  vpt.push_back(0.04250);
    hpt.push_back(0.00700);  vpt.push_back(0.04250);
    hpt.push_back(0.00600);  vpt.push_back(0.04250);
    hpt.push_back(0.00500);  vpt.push_back(0.04250);
    hpt.push_back(0.00400);  vpt.push_back(0.04250);
    hpt.push_back(0.00300);  vpt.push_back(0.04240);
    hpt.push_back(0.00200);  vpt.push_back(0.04230);
    hpt.push_back(0.00100);  vpt.push_back(0.04220);
    hpt.push_back(0.00000);  vpt.push_back(0.04210);
    hpt.push_back(-0.00100); vpt.push_back(0.04190);
    hpt.push_back(-0.00200); vpt.push_back(0.04160);
    hpt.push_back(-0.00300); vpt.push_back(0.04130);
    hpt.push_back(-0.00400); vpt.push_back(0.04090);
    hpt.push_back(-0.00500); vpt.push_back(0.04060);
    hpt.push_back(-0.00600); vpt.push_back(0.04020);
    hpt.push_back(-0.00700); vpt.push_back(0.03975);
    hpt.push_back(-0.00800); vpt.push_back(0.03900);
    hpt.push_back(-0.00900); vpt.push_back(0.03820);
    hpt.push_back(-0.01000); vpt.push_back(0.03680);
    hpt.push_back(-0.01100); vpt.push_back(0.03590);
    hpt.push_back(-0.01200); vpt.push_back(0.03500);
    hpt.push_back(-0.01300); vpt.push_back(0.03400);
    hpt.push_back(-0.01400); vpt.push_back(0.03285);
    hpt.push_back(-0.01500); vpt.push_back(0.03170);
    hpt.push_back(-0.01600); vpt.push_back(0.03000);
    hpt.push_back(-0.01700); vpt.push_back(0.02840);
    hpt.push_back(-0.01800); vpt.push_back(0.02705);
    hpt.push_back(-0.01900); vpt.push_back(0.02500);
    hpt.push_back(-0.01950); vpt.push_back(0.02370);
    hpt.push_back(-0.02000); vpt.push_back(0.02150);
    hpt.push_back(-0.02050); vpt.push_back(0.02020);
    hpt.push_back(-0.02100); vpt.push_back(0.01820);
    hpt.push_back(-0.02150); vpt.push_back(0.01760);
    hpt.push_back(-0.02200); vpt.push_back(0.01240);
    hpt.push_back(-0.02250); vpt.push_back(0.00650);
    hpt.push_back(-0.02270); vpt.push_back(0.00350);
    hpt.push_back(-0.02285); vpt.push_back(0.00160);
    hpt.push_back(-0.02300); vpt.push_back(0.00000);
    hpt.push_back(-0.02285); vpt.push_back(-0.00160);
    hpt.push_back(-0.02270); vpt.push_back(-0.00350);
    hpt.push_back(-0.02250); vpt.push_back(-0.00650);
    hpt.push_back(-0.02200); vpt.push_back(-0.01240);
    hpt.push_back(-0.02150); vpt.push_back(-0.01760);
    hpt.push_back(-0.02100); vpt.push_back(-0.01820);
    hpt.push_back(-0.02050); vpt.push_back(-0.02020);
    hpt.push_back(-0.02000); vpt.push_back(-0.02150);
    hpt.push_back(-0.01950); vpt.push_back(-0.02370);
    hpt.push_back(-0.01900); vpt.push_back(-0.02450);
    hpt.push_back(-0.01800); vpt.push_back(-0.02560);
    hpt.push_back(-0.01700); vpt.push_back(-0.02680);
    hpt.push_back(-0.01600); vpt.push_back(-0.02810);
    hpt.push_back(-0.01500); vpt.push_back(-0.02930);
    hpt.push_back(-0.01400); vpt.push_back(-0.03010);
    hpt.push_back(-0.01300); vpt.push_back(-0.03220);
    hpt.push_back(-0.01200); vpt.push_back(-0.03340);
    hpt.push_back(-0.01100); vpt.push_back(-0.03460);
    hpt.push_back(-0.01000); vpt.push_back(-0.03550);
    hpt.push_back(-0.00900); vpt.push_back(-0.03630);
    hpt.push_back(-0.00800); vpt.push_back(-0.03700);
    hpt.push_back(-0.00700); vpt.push_back(-0.03780);
    hpt.push_back(-0.00600); vpt.push_back(-0.03880);
    hpt.push_back(-0.00500); vpt.push_back(-0.03990);
    hpt.push_back(-0.00400); vpt.push_back(-0.04100);
    hpt.push_back(-0.00300); vpt.push_back(-0.04140);
    hpt.push_back(-0.00200); vpt.push_back(-0.04170);
    hpt.push_back(-0.00100); vpt.push_back(-0.04210);
    hpt.push_back(0.00000);  vpt.push_back(-0.04260);
    hpt.push_back(0.00100);  vpt.push_back(-0.04300);
    hpt.push_back(0.00200);  vpt.push_back(-0.04350);
    hpt.push_back(0.00300);  vpt.push_back(-0.04370);
    hpt.push_back(0.00400);  vpt.push_back(-0.04400);
    hpt.push_back(0.00500);  vpt.push_back(-0.04420);
    hpt.push_back(0.00600);  vpt.push_back(-0.04440);
    hpt.push_back(0.00700);  vpt.push_back(-0.04460);
    hpt.push_back(0.00800);  vpt.push_back(-0.04480);
    hpt.push_back(0.00900);  vpt.push_back(-0.04490);
    hpt.push_back(0.01000);  vpt.push_back(-0.04500);
    hpt.push_back(0.01100);  vpt.push_back(-0.04500);
    hpt.push_back(0.01200);  vpt.push_back(-0.04505);
    hpt.push_back(0.01300);  vpt.push_back(-0.04500);
    hpt.push_back(0.012995); vpt.push_back(-0.04380);
    hpt.push_back(0.01299); vpt.push_back(-0.04220);
    hpt.push_back(0.012985); vpt.push_back(-0.04050);
    hpt.push_back(0.01298); vpt.push_back(-0.03890);
    hpt.push_back(0.012975); vpt.push_back(-0.03720);
    hpt.push_back(0.01297); vpt.push_back(-0.03560);
    hpt.push_back(0.012965); vpt.push_back(-0.03390);
    hpt.push_back(0.01296); vpt.push_back(-0.03230);
    hpt.push_back(0.012955); vpt.push_back(-0.03060);
    hpt.push_back(0.01295); vpt.push_back(-0.02900);
    hpt.push_back(0.012945); vpt.push_back(-0.02740);
    hpt.push_back(0.01294); vpt.push_back(-0.02570);
    hpt.push_back(0.012935); vpt.push_back(-0.02400);
    hpt.push_back(0.01293); vpt.push_back(-0.02240);
    hpt.push_back(0.012925); vpt.push_back(-0.02070);
    hpt.push_back(0.01292); vpt.push_back(-0.01910);
    hpt.push_back(0.012915); vpt.push_back(-0.01740);
    hpt.push_back(0.01291); vpt.push_back(-0.01580);
    hpt.push_back(0.012905); vpt.push_back(-0.01420);
    hpt.push_back(0.01290); vpt.push_back(-0.01250);
    hpt.push_back(0.012895); vpt.push_back(-0.01090);
    hpt.push_back(0.01289); vpt.push_back(-0.00920);
    hpt.push_back(0.01265); vpt.push_back(-0.00808);
    hpt.push_back(0.01230); vpt.push_back(-0.00646);
    hpt.push_back(0.01195); vpt.push_back(-0.00484);
    hpt.push_back(0.01160); vpt.push_back(-0.00324);
    hpt.push_back(0.01140); vpt.push_back(-0.00161);
    hpt.push_back(0.01130); vpt.push_back(-0.00000);
    hpt.push_back(0.01140); vpt.push_back(-0.00161);
    hpt.push_back(0.01160); vpt.push_back(0.00324);
    hpt.push_back(0.01195); vpt.push_back(0.00484);
    hpt.push_back(0.01230); vpt.push_back(0.00646);
    hpt.push_back(0.01265); vpt.push_back(0.00808);
    hpt.push_back(0.01289); vpt.push_back(0.00920);
    hpt.push_back(0.012895); vpt.push_back(0.01070);
    hpt.push_back(0.01290); vpt.push_back(0.01223);
    hpt.push_back(0.012905); vpt.push_back(0.01370);
    hpt.push_back(0.01291); vpt.push_back(0.01525);
    hpt.push_back(0.012915); vpt.push_back(0.01675);
    hpt.push_back(0.01292); vpt.push_back(0.01826);
    hpt.push_back(0.012925); vpt.push_back(0.01975);
    hpt.push_back(0.01293); vpt.push_back(0.02128);
    hpt.push_back(0.012935); vpt.push_back(0.02278);
    hpt.push_back(0.01294); vpt.push_back(0.02430);
    hpt.push_back(0.012945); vpt.push_back(0.02580);
    hpt.push_back(0.01295); vpt.push_back(0.02733);
    hpt.push_back(0.012955); vpt.push_back(0.02883);
    hpt.push_back(0.01296); vpt.push_back(0.03034);
    hpt.push_back(0.012965); vpt.push_back(0.03184);
    hpt.push_back(0.01297); vpt.push_back(0.03337);
    hpt.push_back(0.012975); vpt.push_back(0.034870);
    hpt.push_back(0.01298); vpt.push_back(0.03638);
    hpt.push_back(0.012985); vpt.push_back(0.037880);
    hpt.push_back(0.01299); vpt.push_back(0.03942);
    hpt.push_back(0.012995); vpt.push_back(0.040920);
*/

    // A first try at PREX-2 collimation.  Oct 25, 2020
    //based on PREX-2 new database
    hpt.push_back(0.01500);  vpt.push_back(0.03890);
    hpt.push_back(0.01400);  vpt.push_back(0.03890);
    hpt.push_back(0.01300);  vpt.push_back(0.03890);
    hpt.push_back(0.01200);  vpt.push_back(0.03890);
    hpt.push_back(0.01100);  vpt.push_back(0.03890);
    hpt.push_back(0.01000);  vpt.push_back(0.03890);
    hpt.push_back(0.00900);  vpt.push_back(0.03890);
    hpt.push_back(0.00800);  vpt.push_back(0.03890);
    hpt.push_back(0.00700);  vpt.push_back(0.03890);
    hpt.push_back(0.00600);  vpt.push_back(0.03890);
    hpt.push_back(0.00500);  vpt.push_back(0.03890);
    hpt.push_back(0.00400);  vpt.push_back(0.03890);
    hpt.push_back(0.00300);  vpt.push_back(0.03880);
    hpt.push_back(0.00200);  vpt.push_back(0.03870);
    hpt.push_back(0.00100);  vpt.push_back(0.03860);
    hpt.push_back(0.00000);  vpt.push_back(0.03850);
    hpt.push_back(-0.00100); vpt.push_back(0.03840);
    hpt.push_back(-0.00200); vpt.push_back(0.03820);
    hpt.push_back(-0.00300); vpt.push_back(0.03810);
    hpt.push_back(-0.00400); vpt.push_back(0.03780);
    hpt.push_back(-0.00500); vpt.push_back(0.03760);
    hpt.push_back(-0.00600); vpt.push_back(0.03730);
    hpt.push_back(-0.00700); vpt.push_back(0.03690);
    hpt.push_back(-0.00800); vpt.push_back(0.03640);
    hpt.push_back(-0.00900); vpt.push_back(0.03590);
    hpt.push_back(-0.01000); vpt.push_back(0.03520);
    hpt.push_back(-0.01100); vpt.push_back(0.03430);
    hpt.push_back(-0.01200); vpt.push_back(0.03360);
    hpt.push_back(-0.01300); vpt.push_back(0.03270);
    hpt.push_back(-0.01400); vpt.push_back(0.03160);
    hpt.push_back(-0.01500); vpt.push_back(0.03000);
    hpt.push_back(-0.01600); vpt.push_back(0.02800);
    hpt.push_back(-0.01700); vpt.push_back(0.02760);
    hpt.push_back(-0.01800); vpt.push_back(0.02650);
    hpt.push_back(-0.01900); vpt.push_back(0.02500);
    hpt.push_back(-0.01950); vpt.push_back(0.02370);
    hpt.push_back(-0.02000); vpt.push_back(0.02150);
    hpt.push_back(-0.02050); vpt.push_back(0.02020);
    hpt.push_back(-0.02100); vpt.push_back(0.01820);
    hpt.push_back(-0.02150); vpt.push_back(0.01760);
    hpt.push_back(-0.02200); vpt.push_back(0.01240);
    hpt.push_back(-0.02250); vpt.push_back(0.00650);
    hpt.push_back(-0.02270); vpt.push_back(0.00350);
    hpt.push_back(-0.02285); vpt.push_back(0.00160);
    hpt.push_back(-0.02300); vpt.push_back(0.00000);
    hpt.push_back(-0.02295); vpt.push_back(-0.00160);
    hpt.push_back(-0.02290); vpt.push_back(-0.00350);
    hpt.push_back(-0.02280); vpt.push_back(-0.00650);
    hpt.push_back(-0.02265); vpt.push_back(-0.01240);
    hpt.push_back(-0.02240); vpt.push_back(-0.01900);
    hpt.push_back(-0.02200); vpt.push_back(-0.02050);
    hpt.push_back(-0.02140); vpt.push_back(-0.02200);
    hpt.push_back(-0.02050); vpt.push_back(-0.02600);
    hpt.push_back(-0.01970); vpt.push_back(-0.02950);
    hpt.push_back(-0.01900); vpt.push_back(-0.03060);
    hpt.push_back(-0.01800); vpt.push_back(-0.03200);
    hpt.push_back(-0.01700); vpt.push_back(-0.03350);
    hpt.push_back(-0.01600); vpt.push_back(-0.03450);
    hpt.push_back(-0.01500); vpt.push_back(-0.03510);
    hpt.push_back(-0.01400); vpt.push_back(-0.03600);
    hpt.push_back(-0.01300); vpt.push_back(-0.03700);
    hpt.push_back(-0.01200); vpt.push_back(-0.03800);
    hpt.push_back(-0.01100); vpt.push_back(-0.03900);
    hpt.push_back(-0.01000); vpt.push_back(-0.04000);
    hpt.push_back(-0.00900); vpt.push_back(-0.04090);
    hpt.push_back(-0.00800); vpt.push_back(-0.04170);
    hpt.push_back(-0.00700); vpt.push_back(-0.04240);
    hpt.push_back(-0.00600); vpt.push_back(-0.04300);
    hpt.push_back(-0.00500); vpt.push_back(-0.04360);
    hpt.push_back(-0.00400); vpt.push_back(-0.04410);
    hpt.push_back(-0.00300); vpt.push_back(-0.04450);
    hpt.push_back(-0.00200); vpt.push_back(-0.04500);
    hpt.push_back(-0.00100); vpt.push_back(-0.04570);
    hpt.push_back(0.00000);  vpt.push_back(-0.04620);
    hpt.push_back(0.00100);  vpt.push_back(-0.04650);
    hpt.push_back(0.00200);  vpt.push_back(-0.04680);
    hpt.push_back(0.00300);  vpt.push_back(-0.04700);
    hpt.push_back(0.00400);  vpt.push_back(-0.04740);
    hpt.push_back(0.00500);  vpt.push_back(-0.04780);
    hpt.push_back(0.00600);  vpt.push_back(-0.04805);
    hpt.push_back(0.00700);  vpt.push_back(-0.04840);
    hpt.push_back(0.00800);  vpt.push_back(-0.04870);
    hpt.push_back(0.00900);  vpt.push_back(-0.04900);
    hpt.push_back(0.01000);  vpt.push_back(-0.04925);
    hpt.push_back(0.01100);  vpt.push_back(-0.04950);
    hpt.push_back(0.01200);  vpt.push_back(-0.04970);
    hpt.push_back(0.01300);  vpt.push_back(-0.04980);
    hpt.push_back(0.01400);  vpt.push_back(-0.05000);
    hpt.push_back(0.01500);  vpt.push_back(-0.05020);
    hpt.push_back(0.01498);  vpt.push_back(-0.05000);
    hpt.push_back(0.01496); vpt.push_back(-0.04900);
    hpt.push_back(0.01495); vpt.push_back(-0.04750);
    hpt.push_back(0.01493); vpt.push_back(-0.04600);
    hpt.push_back(0.01491); vpt.push_back(-0.04450);
    hpt.push_back(0.01489); vpt.push_back(-0.04300);
    hpt.push_back(0.01487); vpt.push_back(-0.04250);
    hpt.push_back(0.01485); vpt.push_back(-0.04100);
    hpt.push_back(0.01483); vpt.push_back(-0.03950);
    hpt.push_back(0.01481); vpt.push_back(-0.03800);
    hpt.push_back(0.01479); vpt.push_back(-0.03650);
    hpt.push_back(0.01480); vpt.push_back(-0.03500);
    hpt.push_back(0.01470); vpt.push_back(-0.03350);
    hpt.push_back(0.01460); vpt.push_back(-0.03200);
    hpt.push_back(0.01450); vpt.push_back(-0.03050);
    hpt.push_back(0.01440); vpt.push_back(-0.02880);
    hpt.push_back(0.01430); vpt.push_back(-0.02700);
    hpt.push_back(0.01420); vpt.push_back(-0.02500);
    hpt.push_back(0.01410); vpt.push_back(-0.02300);
    hpt.push_back(0.01405); vpt.push_back(-0.02050);
    hpt.push_back(0.01400); vpt.push_back(-0.01900);
    hpt.push_back(0.01395); vpt.push_back(-0.01750);
    hpt.push_back(0.01390); vpt.push_back(-0.01600);
    hpt.push_back(0.01385); vpt.push_back(-0.01550);
    hpt.push_back(0.01380); vpt.push_back(-0.01420);
    hpt.push_back(0.01370); vpt.push_back(-0.01250);
    hpt.push_back(0.01350); vpt.push_back(-0.01090);
    hpt.push_back(0.01340); vpt.push_back(-0.00920);
    hpt.push_back(0.01330); vpt.push_back(-0.00808);
    hpt.push_back(0.01310); vpt.push_back(-0.00646);
    hpt.push_back(0.01300); vpt.push_back(-0.00484);
    hpt.push_back(0.01300); vpt.push_back(-0.00324);
    hpt.push_back(0.01300); vpt.push_back(-0.00161);
    hpt.push_back(0.01300); vpt.push_back(-0.00000);
    hpt.push_back(0.01300); vpt.push_back(0.00161);
    hpt.push_back(0.01300); vpt.push_back(0.00324);
    hpt.push_back(0.01310); vpt.push_back(0.00484);
    hpt.push_back(0.01330); vpt.push_back(0.00646);
    hpt.push_back(0.01340); vpt.push_back(0.00808);
    hpt.push_back(0.01350); vpt.push_back(0.00920);
    hpt.push_back(0.01370); vpt.push_back(0.01070);
    hpt.push_back(0.01380); vpt.push_back(0.01223);
    hpt.push_back(0.01390); vpt.push_back(0.01370);
    hpt.push_back(0.01400); vpt.push_back(0.01525);
    hpt.push_back(0.01410); vpt.push_back(0.01675);
    hpt.push_back(0.01420); vpt.push_back(0.01826);
    hpt.push_back(0.01430); vpt.push_back(0.01975);
    hpt.push_back(0.01440); vpt.push_back(0.02128);
    hpt.push_back(0.01445); vpt.push_back(0.02278);
    hpt.push_back(0.01450); vpt.push_back(0.02430);
    hpt.push_back(0.01455); vpt.push_back(0.02580);
    hpt.push_back(0.01460); vpt.push_back(0.02733);
    hpt.push_back(0.01465); vpt.push_back(0.02883);
    hpt.push_back(0.01470); vpt.push_back(0.03034);
    hpt.push_back(0.01475); vpt.push_back(0.03184);
    hpt.push_back(0.01480); vpt.push_back(0.03337);
    hpt.push_back(0.01485); vpt.push_back(0.034870);
    hpt.push_back(0.01495); vpt.push_back(0.03638);
    hpt.push_back(0.01496); vpt.push_back(0.03700);
    hpt.push_back(0.01497); vpt.push_back(0.03760);
    hpt.push_back(0.01498); vpt.push_back(0.03810);

    npts = vpt.size();
  
    Float_t osign = 1;  // -1 for L-HRS, +1 for R-HRS (this changed in 2020)
    ohmin = 99999;
    ohmax = -99999;
    ovmin = 99999;
    ovmax = -99999;

    Float_t hshift = 0.0;

    for (Int_t i = 0; i < npts; i++) {
  // Possible Sign flip for an HRS (compared to hamc)
        hpt[i] = osign*hpt[i] + hshift;
        hpt[i] = blowup * hpt[i];
        vpt[i] = blowup * vpt[i];
        if (hpt[i] < ohmin) ohmin = hpt[i];
        if (hpt[i] > ohmax) ohmax = hpt[i];
        if (vpt[i] < ovmin) ovmin = vpt[i];
        if (vpt[i] > ovmax) ovmax = vpt[i];

	//	std::cout << "Emp. angle "<<i<<"  "<<hpt[i]<<"  "<<vpt[i]<<"  signs "<<hsign[i]<<"   "<<vsign[i]<<std::endl;
    }

     Print();
  }
  void Print() const { 
    std::cout<<"Empirical angle collimation. **** Num points = "<<npts<<std::endl;
    for (Int_t i = 0; i < npts; i++) {
      std::cout << "Pt "<<i<<"  "<<hpt[i]<<"  "<<vpt[i];
     }
    hamcAperture::Print();
  }
  Int_t WhichBox(Float_t vert, Float_t horiz) const {
 // Check angles from target.  Units are radians.
 // Here, horiz is horizontal angle, vert = vertical

    if (debug) std::cout << "Into WhichBox vert,horiz =  "<<vert<<"  "<<horiz<<std::endl;
    if (horiz < ohmin || horiz > ohmax) return -1;
    if (vert < ovmin || vert > ovmax) return -1;

    Int_t i, j, c = 0;
    for (i = 0, j = npts-1; i < npts; j = i++) {
      if ( ((vpt[i]>vert) != (vpt[j]>vert)) &&
	   (horiz < (hpt[j]-hpt[i]) * (vert-vpt[i]) / (vpt[j]-vpt[i]) + hpt[i]) )
	c = !c;
    }
    if (c == 0) return -1;
    return 1;

  }
  Bool_t CheckAccept(hamcTrack *trk) const {
    if (!trk) return kFALSE;
    Float_t theta = trk->tvect->GetTheta();
    Float_t phi = trk->tvect->GetPhi();
    Int_t idx = WhichBox(theta,phi);
    if (idx == 1) return kTRUE;
    return kFALSE;
  }
private:
  Int_t npts;
  Int_t debug;
  Float_t blowup, infty;
  Float_t ohmin, ohmax; 
  Float_t ovmin, ovmax; 
  std::vector<Float_t> vpt,hpt;

#ifndef NODICT
ClassDef (hamcAngleCollim, 0)   // Empirical angle collimation.
#endif
};



#endif


   
