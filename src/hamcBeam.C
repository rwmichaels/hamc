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
  beam_current = 100; // uA
  polarization = 0.85;
  polerr = 0.018;  // fractional error (0.01 = 1%).
}

hamcBeam::hamcBeam(Float_t ebeam, Float_t esigma) : hamcTrack("electron"),rastered(kTRUE),xrast(0.4),yrast(0.4)
{
  trktype = "beam";
  E0 = ebeam;
  E0sigma = esigma;
  beam_current = 100; // uA
  polarization = 0.85;
  polerr = 0.018;  // fractional error (0.01 = 1%).
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
     if (parser.IsFound("current")) {
       beam_current = parser.GetData();
     }   
     if (parser.IsFound("polarization")) {
       polarization = parser.GetData();
     }   

     cout << "Beam values   E0  = "<<E0<<"   E0sigma "<<E0sigma<<endl;
     if (IsRastered()) {
       cout << "Raster  x "<<xrast<<"  y "<<yrast<<endl;
     }

     cout << "Beam current "<<beam_current<<" uA "<<endl;

     dx_iter = 0;
     dy_iter = 0;
     dE_iter = 0;
     dtheta_iter = 0;
     dphi_iter = 0;

     parser.Load(expt->inout->GetStrVect("kick:track"));
     parser.Print();
     if (parser.IsFound("x")){
       if (expt->inout->numiter > 1) dx_iter = parser.GetData();
       cout << "Will iterate x by "<<dx_iter<<" meters"<<endl;
     }
     if (parser.IsFound("y")){
       if (expt->inout->numiter > 1) dy_iter = parser.GetData();
       cout << "Will iterate y by "<<dy_iter<<" meters"<<endl;
     }
     if (parser.IsFound("E")){
       if (expt->inout->numiter > 1) dE_iter = parser.GetData();
       cout << "Will iterate E by "<<dE_iter<<"  GeV "<<endl;
     }
     if (parser.IsFound("theta")){
       if (expt->inout->numiter > 1) dtheta_iter = parser.GetData();
       cout<< "Will iterate the beam incident angle theta by "<<dtheta_iter<<" radians"<<endl;
     }
     if (parser.IsFound("phi")){
       if (expt->inout->numiter > 1) dphi_iter = parser.GetData();
       cout<< "Will iterate the beam incident angle phi by "<<dphi_iter<<" radians"<<endl;
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

  Float_t x_position, y_position;

  tvect->Clear();


  if (IsRastered()) { 
     x_position = -0.5*xrast + xrast*gRandom->Rndm(1);
     y_position = -0.5*yrast + yrast*gRandom->Rndm(1);

    if (expt->iteration == 1){
      x_position += dx_iter;
      y_position += dy_iter;
    }

    tvect->PutX(x_position);
    tvect->PutY(y_position);


  } else {
    tvect->Clear();
  }

  tvect->PutZ(expt->target->GetZScatt());

  energy = E0 + E0sigma*gRandom->Gaus();
  if (expt->iteration == 1)
    energy += dE_iter;

  Radiate(expt);  // modifies the energy

  if (energy < mass) energy = mass;  // extrema of rad tail

  pmom = TMath::Sqrt(energy*energy - mass*mass);

  if (expt->iteration == 1) {

    if (dtheta_iter > 1.E-8){
    
      pnoms_x = pmom*TMath::Sin(dtheta_iter);
      pnoms_y = 0;
      pnoms_z = pmom*TMath::Cos(dtheta_iter);
  
   } else if (dphi_iter > 1.E-8) {
    
      pnoms_x = 0;
      pnoms_y = pmom*TMath::Sin(dphi_iter);
      pnoms_z = pmom*TMath::Cos(dphi_iter);
    }
  
  } else {

    // Assume the beam is along the Z axis.
    
    pnoms_x = 0;
    pnoms_y = 0;
    pnoms_z = pmom;
  }

 // Apply multiple scattering

  /*
   * spr:  antiquated
  Float_t radlen;
  Float_t tlen =  expt->target->GetLength();  // "L" of target
// zscat varies from -L/2 to +L/2
  Float_t zscat = expt->target->GetZScatt();
  Float_t x0rad = expt->target->GetRadLength();
  radlen = 0;
  if (x0rad != 0) radlen = (((tlen/2)+zscat)/tlen) * x0rad;
  MultScatt(radlen, ITARGET);
  */

  //  cout << "\n ******* Calling Mult Scatt in beam "<<endl;
   MultScatt(expt->target->GetPartialScatt(pmom), ITARGET);    // this modifies/updates plab (momentum)
                                 // Note, MultScatt is a member of hamcTrack
                                 // and hamcBeam inherits from hamcTrack.

  //  cout << "Beam Pmom "<<endl;
  //  cout << "No MS "<<pnoms_x<<"   "<<pnoms_y<<"   "<<pnoms_z<<endl;
  //  cout << "With MS "<<plab_x<<"   "<<plab_y<<"   "<<plab_z<<endl;


  return OK;
}

Int_t hamcBeam::Radiate(hamcExpt* expt) {

// The internal is for 1/2 t_equivalent, so it's
// what to subtract before and after scattering.

  Float_t dE = expt->physics->eloss->GetDeExternIn() + 
               expt->physics->eloss->GetDeIntern() +
               expt->physics->eloss->GetDeIonIn();

  energy = energy - dE;

  return 0;

}





   
