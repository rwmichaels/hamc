#ifndef ROOT_hamcTransMat
#define ROOT_hamcTransMat

//  hamcTransMat -- 1st order Transport matrix for HRS
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcTrans.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <string>
#include <map>

class hamcSpecHRS;

// Utility class
class hamcQuad {
// Quadrupole beamline element.
//         (x,x')    cm/rad  (not mrad)
//         (x',x)    rad/cm
 public:
  hamcQuad(Float_t bbp, Float_t rrad, Float_t llen) : 
// bbp = field (gauss) at pole,  rrad = radius of quad (cm), llen = length (cm)
    bp(bbp), rad(rrad), len(llen), focusx(1), pcut(0.01), huge(999999), debug(0) 
    {  grad = bp / rad;
      kappa0 = TMath::Sqrt(4.8e-10 * grad / 1.6e-3); }
  virtual ~hamcQuad() {};
  void SetFocusX() { focusx = 1; };
  void SetFocusY() { focusx = 0; };
  Float_t GetMatrix(Float_t pmom, Int_t index) {
    if (pmom < pcut) return huge;
    kappa = kappa0 / TMath::Sqrt(pmom);
    kL = kappa * len;
    exp1 = exp(kL);
    exp2 = exp(-1*kL);
    coshx = (exp1 + exp2)/2;
    sinhx = (exp1 - exp2)/2;
     if (debug) 
       std::cout << "focus "<<focusx<<"  pmom "<<pmom<<"  kappa "<<kappa<<"  kL "<<kL<<std::endl;
    if (focusx) {
      if (index == 0) return TMath::Cos(kL);
      if (index == 1) return (1/kappa)*TMath::Sin(kL);
      if (index == 2) return -1*kappa*TMath::Sin(kL);
      if (index == 3) return TMath::Cos(kL);
      if (index == 4) return coshx;
      if (index == 5) return (1/kappa)*sinhx;
      if (index == 6) return kappa*sinhx;
      if (index == 7) return coshx;
    } else {
      if (index == 0) return coshx;
      if (index == 1) return (1/kappa)*sinhx;
      if (index == 2) return kappa*sinhx;
      if (index == 3) return coshx;
      if (index == 4) return TMath::Cos(kL);
      if (index == 5) return (1/kappa)*TMath::Sin(kL);
      if (index == 6) return -1*kappa*TMath::Sin(kL);
      if (index == 7) return TMath::Cos(kL);
    }
    return 0;
  }
  Float_t GetLen() { return len; };
  void Print() {
    std::cout << "Quad values "<<std::endl;
    std::cout << "grad "<<grad<<"    kappa "<<kappa<<std::endl;
    if (focusx) {
      std::cout<<"focus in X"<<std::endl;
    } else {
      std::cout<<"focus in Y"<<std::endl;
    }
  }
 private:
  Float_t bp, rad, len;
  Int_t focusx;
  Float_t pcut, huge;
  Int_t debug;
  Float_t grad, kappa0;
  Float_t exp1, exp2;
  Float_t kappa, kL, coshx, sinhx;
#ifndef NODICT
ClassDef (hamcQuad, 0)   // Quadrupole beamline element
#endif
};


class hamcMatrix {  // 5x5 matrix for transport
 public: 
  hamcMatrix(Float_t x00, Float_t x01, Float_t x02, Float_t x03, Float_t x04, 
             Float_t x10, Float_t x11, Float_t x12, Float_t x13, Float_t x14, 
             Float_t x20, Float_t x21, Float_t x22, Float_t x23, Float_t x24, 
             Float_t x30, Float_t x31, Float_t x32, Float_t x33, Float_t x34, 
             Float_t x40, Float_t x41, Float_t x42, Float_t x43, Float_t x44) {
    matrix = new Float_t[25];
    matrix[0]  = x00; matrix[1]  = x01; matrix[2]  = x02; matrix[3]  = x03; matrix[4]  = x04;
    matrix[5]  = x10; matrix[6]  = x11; matrix[7]  = x12; matrix[8]  = x13; matrix[9]  = x14;
    matrix[10] = x20; matrix[11] = x21; matrix[12] = x22; matrix[13] = x23; matrix[14] = x24;
    matrix[15] = x30; matrix[16] = x31; matrix[17] = x32; matrix[18] = x33; matrix[19] = x34;
    matrix[20] = x40; matrix[21] = x41; matrix[22] = x42; matrix[23] = x43; matrix[24] = x44;
  }
  Float_t Get(Int_t row, Int_t col) {  // row = 0,1,2,3,4    col = 0,1,2,3,4
    Int_t index = 5*row+col;
    if (index < 0 || index > 24) return 0;
    return matrix[index];
  }
  void Print() {
    std::cout << "Matrix :"<<std::endl;
    for (Int_t row=0; row<5; row++) {
      for (Int_t col=0; col<5; col++) {
	std::cout << Get(row,col) << "  ";
        if (col==4) std::cout<<std::endl;
      }
    }
    std::cout<<std::endl;
  }
  virtual ~hamcMatrix() { delete matrix;};
private:
   Float_t *matrix;
#ifndef NODICT
ClassDef (hamcMatrix, 0)   // One Transport Matrix
#endif
};
  
 
class hamcTransMat : public hamcTrans {

  public:

     hamcTransMat();
     virtual ~hamcTransMat(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

  private: 

     Int_t MatrixMult(Int_t idx, hamcTransVect& tvec) const;  
     std::vector<hamcMatrix *> matrix;
     std::vector<Int_t> ptr_matrix;

     hamcTransMat(const hamcTransMat& trans);
     hamcTransMat& operator=(const hamcTransMat& trans);

#ifndef NODICT
ClassDef (hamcTransMat, 0)   // Matrix Model for HRS Transport
#endif
};

#endif



   
