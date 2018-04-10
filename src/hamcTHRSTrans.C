//  hamcTHRSTrans -- Wrapper class for the THRSTrans library
//  which is the separately maintained software Seamus Riordan wrote.
//  R. Michaels, Dec 2017

#include "hamcTHRSTrans.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTHRSTrans)
#endif


hamcTHRSTrans::hamcTHRSTrans() 
{
  trans=0;
}

hamcTHRSTrans::~hamcTHRSTrans() {
}

Int_t hamcTHRSTrans::Init(hamcSpecHRS *spect) {

  which_spectrom = spect->which_spectrom;
  collim_distance = spect->GetCollimDist();

  double q1,q2,q3,k1,k2;
  q1 = spect->quad1;
  q2 = spect->quad2;
  q3 = spect->quad3;
  k1 = spect->dipk1;
  k2 = spect->dipk2;
  fTune = spect->fTune;

  trans = new THRSTrans(q1, q2, q3, k1, k2, fTune);

  cout << "\nUsing THRSTrans  for tune "<<DoInterpret(fTune) << endl << endl;
  cout << "Quads "<<q1<<"  "<<q2<<"  "<<q3<<"   dipole "<<k1<<"  "<<k2<<endl;
  trans->ShowOutput();

  return OK;

}

void hamcTHRSTrans::Print() {
}

TString hamcTHRSTrans::DoInterpret( THRSTrans::tune_t thetune) {
  TString result=" Undefined ";
  if(thetune==THRSTrans::kStd)  result = " Standard HRS ";
  if(thetune==THRSTrans::kPREX) result = " PREX ";
  if(thetune==THRSTrans::kCREX) result = " CREX ";
  if(thetune==THRSTrans::kAPEX) result = " APEX (??)  ";
  return result;
}

Int_t hamcTHRSTrans::TransForm(hamcTrack *trk, Int_t where) const {

// Transforms a track "trk" to the destination = where.
// The origin of the track is known (trk->origin = ITARGET, ICOLLIM, etc).
// Seamus's class THRSTrans is used as an underlying library.
//
// Input and output is transport vector. +X is down.
// The (x,y) is not correct for the aperture functions from LeRose, though;
// aperture functions in Transport coordinates should be used.

// NOTE: For warm septum there is no distinction between Left
// and Right HRS at the moment.

  if (!trans) {
    cout << "ERROR: hamcTHRSTrans: No THRSTrans object ! \n"<<endl;
    return 0;
  }

  TVectorD iv(20),v(20);
  TMatrixD hrsmatrix(20,20);
  Int_t idx=0;

  Int_t debug_trans = 0; // to get do some debugging here.
  Int_t origin = trk->origin;

  iv[THRSTrans::kX]  = trk->tvect_orig->GetX();
  iv[THRSTrans::kTh] = trk->tvect_orig->GetTheta();
  iv[THRSTrans::kY]  = trk->tvect_orig->GetY();
  iv[THRSTrans::kPh] = trk->tvect_orig->GetPhi();
  iv[THRSTrans::kd]  = trk->tvect_orig->GetDpp();

  trans->fillvector(iv);  // add 2nd order

  if (origin == ITARGET) {

    switch(where) {

      case ISEPTIN: 

        hrsmatrix = *(trans->GetTransport(2));
        break;

      case ISEPTOUT: 
 
        hrsmatrix = *(trans->GetTransport(4));
        break;

      case ICOLLIM: 
      case ICOLLIM2: 
      case ICOLLIM3: 

        hrsmatrix = *(trans->GetTransport(5));
        break;

      case IQ1IN: 

        idx=6;
        if (fTune==THRSTrans::kStd) idx=1;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IQ1EXIT: 

        idx=7;
        if (fTune==THRSTrans::kStd) idx=2;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IQ2IN: 

        idx=8;
        if (fTune==THRSTrans::kStd) idx=3;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IQ2EXIT: 

        idx=9;
        if (fTune==THRSTrans::kStd) idx=4;

        hrsmatrix = *(trans->GetTransport(9));
        break;

      case IDIPIN:

        idx=10;
        if (fTune==THRSTrans::kStd) idx=5;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IDIPEXIT:

        idx=11;
        if (fTune==THRSTrans::kStd) idx=6;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IQ3IN:

        idx=12;
        if (fTune==THRSTrans::kStd) idx=7;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IQ3EXIT:

        idx=13;
        if (fTune==THRSTrans::kStd) idx=8;

        hrsmatrix = *(trans->GetTransport(idx));
        break;

      case IFOCAL:

        hrsmatrix = *(trans->GetTransport());
        break;

      default:

        cout << "hamcTHRSTrans::WARNING: undefined destination "<<where<<endl;
        return 0;

    }

  } else {

        cout << "hamcTHRSTrans::WARNING: undefined origin "<<origin<<endl;
        return 0;
  
  }

  v = (hrsmatrix)*iv;

  trk->tvect->PutX(v[THRSTrans::kX]);
  trk->tvect->PutTheta(v[THRSTrans::kTh]);
  trk->tvect->PutY(v[THRSTrans::kY]);
  trk->tvect->PutPhi(v[THRSTrans::kPh]);
  trk->tvect->PutDpp(v[THRSTrans::kd]);
  
  return OK;

}

