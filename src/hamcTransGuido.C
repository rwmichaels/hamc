//  hamcTransGuido   -- Guido's Model for HRS Transport
//  R. Michaels  Mar 2009

#include "hamcTransGuido.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTransGuido)
#endif


hamcTransGuido::hamcTransGuido() : hamcTrans() 
{
}

hamcTransGuido::~hamcTransGuido() {

}

Int_t hamcTransGuido::Init(hamcSpecHRS *spect) {

// Initialize the parameters

// Notation:  "xDth" means how "x" depends on "th".
// If you have "xDthth" it is second order (multiplies th*th).

// EXAMPLE:

  xDth = 12.5;     // cm/mr
  xDthth = 1.4;
  xDph = 0.4;
  xDphph = 0.66;
  yDth = 0.65;
  yDph = 0.89;
  yDphph = 0.0124;
  thDx = 0.554;
  phDx = 0.343;
  phDxx = 0.044;
  phDy = 0.8;

  return OK;

}

void hamcTransGuido::Print() {
  cout << "Parameters of Guido's transport"<<endl;
  // put print of parameters here
}

Int_t hamcTransGuido::TransForm(hamcTrack *trk, Int_t where) const {
// Transforms "trk" to the focal plane.

  hamcTransVect tvect;

  tvect.PutX(trk->tvect_orig->GetX());
  tvect.PutTheta(trk->tvect_orig->GetTheta());
  tvect.PutY(trk->tvect_orig->GetY());
  tvect.PutPhi(trk->tvect_orig->GetPhi());
  tvect.PutDpp(trk->tvect_orig->GetDpp());
		 		 
  TransForm(&tvect, where);
  
  *trk->tvect = tvect;  // Note, this overwrites the main tvect.

  return OK;

}

Int_t hamcTransGuido::TransForm(hamcTransVect *tvect, Int_t where) const {
// This function, rather than the above, is what is actually used
// by hamc at the moment (Mar 2009).  The reason is we want to get 
// an alternative tvect at the focal plane to compare with other
// models. It's also necessary because other models (e.g. LeRose)
// has the required acceptance function, while Guido's does not.

  if (where != IFOCAL) return 0;

  Float_t xorig,yorig,thorig,phorig,dpporig;
  Float_t x,y,th,ph,dpp;

//  Original vector (typically at target)
  xorig   = tvect->GetX();
  yorig   = tvect->GetY();
  thorig  = tvect->GetTheta();
  phorig  = tvect->GetPhi();
  dpporig = tvect->GetDpp();

  // Dustin, you have to make this what you want.

  x = (xDth * thorig)  + (xDthth * thorig*thorig) + 
    (xDph * phorig) + (xDphph * phorig *phorig);

  y = (yDth * thorig) + (yDph * phorig) + (yDphph * phorig*phorig);

  th = 2*thorig; //(thDx * xorig);  

  ph = 4*phorig; //(phDx * xorig) + (phDxx * xorig*xorig) + (phDy * yorig);

  dpp = dpporig;

// Update the vector
  tvect->Load(x,th,y,ph,dpp);

  return OK;

}

