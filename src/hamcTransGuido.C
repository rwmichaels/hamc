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

//Here are Guido's Coeff from run 2373 (March 19, 2009, Dustin McNulty)
//These are meant for LHRS, but are applied to both arms for now
  yDph = 1.11; // m/r
  yDdpp = -0.29; // m
  yDth = -0.05; // m/r
  yDthth = 1.9; // m/r/r
  yDththth = -6.9; // m/r/r/r
  yDphdpp = 8.6; // m/r
  yDphphdpp = 52.6; // m/r/r
  yDy0 = -0.27; // dimensionless

  xDph = 0.18; // m/r
  xDdpp = 13.20; // m
  xDdppdpp = 18.93; // m
  xDdppdppdppdpp = 350.3; // m
  xDth = 0.12; // m/r
  xDththth = 20; // m/r/r/r
  xDthththth = 1300; // m/r/r/r/r
  xDththththth = 9000; // m/r/r/r/r/r
  xDthdpp = -5.4; // m/r
  xDthththph = 3400; // m/r/r/r/r
  xDththththththph = 37060000; // m/r^7
  xDthdppdpp = 175; // m/r
  xDthththdppdpp = 60000; // m/r/r/r
  xDthththdppdppdpp = 1000000; // m/r/r/r
  xDththththdppdppdpp = 6000000; // m/r/r/r/r
  xDphdpp = -7.625; // m/r
  xDphdppdppdpp = 3121; // m/r
  xDthphdpp = 92.21; // m/r/r
  xDthy0y0 = -100; // 1/r/m
  xDthy0y0y0 = 40000; // 1/r/m/m
  xDdppy0 = -3.2; // dimensionless
  xDy0 = 0.1; // dimensionless

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

  //x = (xDth * thorig)  + (xDthth * thorig*thorig) + (xDph * phorig) + (xDphph * phorig *phorig);

  //y = (yDth * thorig) + (yDph * phorig) + (yDphph * phorig*phorig);

  //th = 2*thorig; //(thDx * xorig);  

  //ph = 4*phorig; //(phDx * xorig) + (phDxx * xorig*xorig) + (phDy * yorig);

  //dpp = dpporig;

  //For Guido's Coeff:
  x = (xDph*phorig)+(xDdpp*dpporig)+(xDdppdpp*dpporig*dpporig)+(xDdppdppdppdpp*dpporig*dpporig*dpporig*dpporig)+(xDth*thorig)+(xDththth*thorig*thorig*thorig)+(xDthththth*thorig*thorig*thorig*thorig)+(xDththththth*thorig*thorig*thorig*thorig*thorig)+(xDthdpp*thorig*dpporig)+(xDthththph*thorig*thorig*thorig*phorig)+(xDththththththph*thorig*thorig*thorig*thorig*thorig*thorig*phorig)+(xDthdppdpp*thorig*dpporig*dpporig)+(xDthththdppdpp*thorig*thorig*thorig*dpporig*dpporig)+(xDthththdppdppdpp*thorig*thorig*thorig*dpporig*dpporig*dpporig)+(xDththththdppdppdpp*thorig*thorig*thorig*thorig*dpporig*dpporig*dpporig)+(xDphdpp*phorig*dpporig)+(xDphdppdppdpp*phorig*dpporig*dpporig*dpporig)+(xDthphdpp*thorig*phorig*dpporig)+(xDthy0y0*thorig*yorig*yorig)+(xDthy0y0y0*thorig*yorig*yorig*yorig)+(xDdppy0*dpporig*yorig)+(xDy0*yorig);

  y = (yDph*phorig)+(yDdpp*dpporig)+(yDth*thorig)+(yDthth*thorig*thorig)+(yDththth*thorig*thorig*thorig)+(yDphdpp*phorig*dpporig)+(yDphphdpp*phorig*phorig*dpporig)+(yDy0*yorig);
  y *= -1.0; //Needed this to match data ?Not sure why yet?

  th = thorig;  

  ph = phorig;  

  dpp = dpporig;

// Update the vector
  tvect->Load(x,th,y,ph,dpp);

  return OK;

}

