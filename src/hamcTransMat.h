#ifndef ROOT_hamcTransMat
#define ROOT_hamcTransMat

//  hamcTransMat -- 1st order Transport matrix for HRS
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <string>
#include <map>

class hamcSpecHRS;

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



   
