//  hamcTransMat   -- Transport model using TRANSPORT matrices
//  R. Michaels  June 2008

#include "hamcTransMat.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTransMat)
#endif


hamcTransMat::hamcTransMat() : hamcTrans() 
{
}

hamcTransMat::~hamcTransMat() {
   for (std::vector<hamcMatrix*>::iterator ith = matrix.begin();
     ith != matrix.end(); ith++) delete *ith;
}

Int_t hamcTransMat::Init(hamcSpecHRS *spect) {

// Initialize the "matrix" and "ptr_matrix" vectors.
// ptr_matrix points to the matrix at a location like ICOLLIM, etc.
// Probably want a multimap instead of parallel vectors ... eventually.

// Traditional Units: cm, mrad, (x=vertical), deltap in percent.

  did_init = kTRUE;

// Units of collim_distance must be meters.
  Float_t collim_distance = spect->GetCollimDist();

  cout << "hamcTransMat:: collim_distance = "<<collim_distance<<endl;

// Yeah, this is wrong if collimator is downstream of septum
// (Will fix it later.)
  matrix.push_back(new hamcMatrix( 
      	         1, 0.1*collim_distance, 0, 0, 0, 
                 0, 1, 0, 0, 0, 
                 0, 0, 1, 0.1*collim_distance, 0,  
                 0, 0, 0, 1, 0, 
                 0, 0, 0, 0, 1));

  ptr_matrix.push_back(ICOLLIM);

  matrix.push_back(new hamcMatrix(   // 2nd collimator
      	         1, 0.1*collim_distance, 0, 0, 0, 
                 0, 1, 0, 0, 0, 
                 0, 0, 1, 0.1*collim_distance, 0,  
                 0, 0, 0, 1, 0, 
                 0, 0, 0, 0, 1));

  ptr_matrix.push_back(ICOLLIM2);


  if (spect->IsColdSeptum()) {           // Cold septum (ca 2005)

     matrix.push_back(new hamcMatrix(
        -2.81214,  0.00   , 0.00   , 0.00   ,  14.06069, 
        -3.18504, -0.3556 , 0.00   , 0.00   ,  24.68781,
          0.00   , 0.00   , 1.01244, 0.04074,  0.12918,
          0.00   , 0.00   ,12.80820, 1.50307,  0.51663,
         2.46416, -0.500  , 0.11314, 0.01731, -0.63948 ));

     ptr_matrix.push_back(IFOCAL);

  }  else if (spect->IsWarmSeptum()) {   // The Warm septum, ca 2008

      matrix.push_back(new hamcMatrix(
         -2.507,  -0.0148, -0.089,   0.0,     14.162,
         -0.309,  -0.401,   0.0,     0.002,    2.501,
          0.008,   0.013,   0.321,  -2.181,   -0.351,
          0.007,   0.008,   0.610,  -1.030,   -0.270,
  	  0.0,     0.0,     0.0,     0.0,       0.0)); 

      ptr_matrix.push_back(IFOCAL);

  } else {   // Default:  Just the HRS and no septum

    matrix.push_back(new hamcMatrix(
        -2.07406,-0.03871, 0.00000, 0.00000, 11.91071,  
        -0.82456,-0.49754, 0.00000, 0.00000, 19.61350,
         0.00000, 0.00000,-0.62335,-0.08311,  0.00000,
         0.00000, 0.00000, 3.81145,-1.09605,  0.00000,
         3.08585,-0.51668, 0.00000, 0.00000, -0.63768 ));

     ptr_matrix.push_back(IFOCAL);

  }

  return OK;

}

void hamcTransMat::Print() {
  cout << "Matrix Transport model "<<endl;
  cout << "Number of matrices "<<matrix.size()<<endl;
  for (vector<hamcMatrix *>::iterator im=matrix.begin(); im != matrix.end(); im++) (*im)->Print();
}

Int_t hamcTransMat::TransForm(hamcTrack *trk, Int_t where) const {

// Transforms a track "trk" to where

// First see if the matrix to this location is defined

  Int_t idx = -1;
  for (Int_t i=0; i<(Int_t)ptr_matrix.size(); i++) {
    if (where == ptr_matrix[i]) idx = i;
  }
  if (idx < 0) return OK;        

  *trk->tvect = *trk->tvect_orig;  // will apply to the track origin.

  return MatrixMult(idx, *trk->tvect);

}

Int_t hamcTransMat::MatrixMult(Int_t idx, hamcTransVect& tvect) const {

#ifdef DEBUG1
  cout << "Mult for "<<idx<<endl;
  tvect.Print();
#endif

  if (idx < 0 || idx > (Int_t)matrix.size()-1) return ERROR;

  vector<Float_t> din, dout;

  din = tvect.Get();
  dout.clear();

  for (Int_t row = 0; row < 5; row++) {

    Float_t result=0;

    for (Int_t col = 0; col < 5; col++) {

       result += matrix[idx]->Get(row,col)*din[col];

    }

    dout.push_back(result);

  }

  tvect.Load(dout);  

#ifdef DEBUG1
  cout << "Result "<<endl;
  tvect.Print();
  cout << "------------------------------"<<endl;
#endif

  return OK;

}
