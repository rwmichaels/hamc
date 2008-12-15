//  hamcBeam   -- Incoming beam
//  R. Michaels  Nov 2008

#include "hamcBeam.h"
#include "hamcTrack.h"
#include "hamcExpt.h"
#include "hamcPhysics.h"
#include "hamcEloss.h"
#include "hamcInout.h"
#include "hamcTarget.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcBeam)
#endif


hamcBeam::hamcBeam() : hamcTrack("electron"),rastered(kTRUE),xrast(0.4),yrast(0.4)
{
  trktype = "beam";
}

hamcBeam::hamcBeam(Float_t ebeam, Float_t esigma) : hamcTrack("electron"),rastered(kTRUE),xrast(0.4),yrast(0.4)
{
  trktype = "beam";
  E0 = ebeam;
  E0sigma = esigma;
}

hamcBeam::~hamcBeam() {
}
   

Int_t hamcBeam::Init(hamcExpt *expt) {


  // Default raster info are defined in the constructor.
  // Can over-ride here with setup file.

     hamcStrParser parser;
     parser.Load(expt->inout->GetStrVect("beam"));
     //     parser.Print();
     if (parser.IsFound("rastered")) {
       rastered = kTRUE;
     }
     if (parser.IsFound("xraster")) {
       xrast = parser.GetData(); 
     }   
     if (parser.IsFound("yraster")) {
       yrast = parser.GetData();
     }
     if (parser.IsFound("E0")) {
       E0 = parser.GetData();
     }
     if (parser.IsFound("E0sigma")) {
       E0sigma = parser.GetData();
     }   

     cout << "Beam values   E0  = "<<E0<<"   E0sigma "<<E0sigma<<endl;
     if (IsRastered()) {
       cout << "Raster  x "<<xrast<<"  y "<<yrast<<endl;
     }

     if (E0 > mass) P0 = TMath::Sqrt(E0*E0 - mass*mass);

     did_init = kTRUE;

     return OK;      
}

Int_t hamcBeam::SetRaster(Float_t xmax, Float_t ymax) {
  if (xmax > 0 && ymax > 0) {
     rastered = kTRUE;
  } else {
     rastered = kFALSE;
  }
  xrast = xmax;
  yrast = ymax;
  return OK;
}

void hamcBeam::Print() {
  cout << "\nSetup of "<<trktype<<endl;
  cout << "pid = "<<pid<<endl;
  if (IsRastered()) {
    cout << "Rastered,  x = "<<xrast<<"  y = "<<yrast<<endl;
  }
  cout << "energy "<<E0<<"     spread "<<E0sigma<<endl;
}


Float_t hamcBeam::GetRaster(Int_t ixy) {
  if (ixy==0) return xrast;
  if (ixy==1) return yrast;
  return 0;
}

Int_t hamcBeam::Generate(hamcExpt *expt) {

// assume the angles don't change, only X,Y due to raster.

  if (!did_init) {
    cout << "hamcBeam::ERROR: uninitialized !"<<endl;
    return ERROR;
  }

// xrast,yrast must in meters.  Likewise target's ZScatt.
// Not that it matters much, but X is vertical in
// Transport coordinates.

  if (IsRastered()) {
    tvect->PutX(-0.5*xrast + xrast*gRandom->Rndm(1));
    tvect->PutY(-0.5*yrast + yrast*gRandom->Rndm(1));
  } else {
    tvect->Clear();
  }

  tvect->PutZ(expt->target->GetZScatt());

  energy = E0 + E0sigma*gRandom->Gaus();
  Radiate(expt);  // modifies the energy

  if (energy < mass) energy = mass;  // extrema of rad tail
  pmom = TMath::Sqrt(energy*energy - mass*mass);

// Assume the beam is along the Z axis.

  plab_x = 0;
  plab_y = 0;
  plab_z = pmom;

  return OK;
}

Int_t hamcBeam::Radiate(hamcExpt* expt) {

// The internal is for 1/2 t_equivalent, so it's
// what to subtract before and after scattering.

  Float_t dE = expt->physics->eloss->GetDeExternIn() + 
               expt->physics->eloss->GetDeIntern() +
               expt->physics->eloss->GetDeIonIn();

  energy = energy - dE;

}





   
