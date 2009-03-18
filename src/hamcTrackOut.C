//  hamcTrackOut -- the track going out from target
//  R. Michaels  July 2008

#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcBeam.h"
#include "hamcPhysics.h"
#include "hamcKine.h"
#include "hamcEvent.h"
#include "hamcInout.h"
#include "hamcSpecHRS.h"
#include "hamcExpt.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTrackOut)
#endif


hamcTrackOut::hamcTrackOut() : hamcTrack("electron"),thetamin(0),thetamax(0),phimin(0),phimax(0),which_hrs(0)
{
  did_init = kFALSE;
  det_dist = 1.2; // meters
  trktype = "out";
}


hamcTrackOut::hamcTrackOut(string pid, Float_t ee, Float_t x, Float_t th, Float_t y, Float_t ph, Float_t dp) : hamcTrack(pid,ee,x,th,y,ph,dp)
{
  did_init = kFALSE;
  trktype = "out";
  which_hrs = 0;
}

hamcTrackOut::~hamcTrackOut() {
}
   

Int_t hamcTrackOut::Init(Int_t ispec, hamcExpt *expt) {

  hamcStrParser parser;

  hamcSpecHRS *spec = expt->GetSpectrom(ispec);
  if (!spec) {
    cout << "hamcTrackOut::ERROR: no spectrometer "<<ispec<<endl;
    return ERROR;
  }

  P0 = spec->GetP0();
  P0sigma = spec->GetP0Sigma();
  which_hrs = spec->which_spectrom;  // affects angle convention

  hamcTrack::Init();

// theta_central, converted to radians.
  theta_central = PI * spec->GetScattAngle() / 180;

// Ranges of angles for the simulation

// Defaults if not defined in input file
   
   thetamin = 0.1*theta_central;
   thetamax = 3.0*theta_central;
   phimin = 0;
   phimax = PI;

// Obtain the ranges from the input file if it exists

   parser.Load(expt->inout->GetStrVect("out_angles"));
   if (parser.IsFound("thetamin")) {
    parser.Print();
    thetamin = PI*parser.GetData()/180; 
   }    
  if (parser.IsFound("thetamax")) {
     thetamax = PI*parser.GetData()/180; 
   }
   if (parser.IsFound("phimin")) {
     phimin = PI*parser.GetData()/180; 
   }
   if (parser.IsFound("phimax")) {
     phimax = PI*parser.GetData()/180; 
   }

   cout << "Track-out angles ranges: "<<thetamin<<"  "<<thetamax<<"  "<<phimin<<"  "<<phimax<<endl;

   parser.Load(expt->inout->GetStrVect("detector"));
   if (parser.IsFound("dist")) {
    det_dist = parser.GetData(); 
    cout << "Z dist. of detector to foc. plane = "<<det_dist<<" m"<<endl;
   }    


   Int_t nbin=120;


// Book a few default histograms 
//
//                        to-weight        brk point
//                           ^   care-acc?   ^     hist id
//                           ^      ^        ^       ^
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "th1", 
	        "Theta Generated", &theta, 120,0.8*thetamin,1.2*thetamax);
//                  ^                ^      ^         ^           ^
//                title           variable  nbin     x-low     x-high
//                               (pointer)

   expt->inout->BookHisto(kTRUE, kFALSE, ITARGET, "th2", 
		"Theta at target weighted by crsec", 
                &theta, nbin, 0.8*thetamin,1.2*thetamax);
   expt->inout->BookHisto(kTRUE, kTRUE, ITARGET, "th3", 
		"Theta at target weighted by crsec, in accept", 
                &theta, nbin, 0.8*thetamin,1.2*thetamax);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "th4", 
		"Theta at focal plane weighted by crsec, in accept", 
                &theta, nbin, 0.8*thetamin,1.2*thetamax);

   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "ph", 
		"Phi at target", &phi, nbin, -0.2,1.2*phimax);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "mom", 
	        "Momentum in HRS", &pmom, nbin, 0.5*P0,1.05*P0);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "momu", 
	        "Momentum in HRS", &pmom, nbin, 0.5*P0,1.1*P0);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "momz", 
	   "Momentum in HRS (zoom)", &pmom, nbin, 0.95*P0,1.02*P0);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "momb", 
	   "Momentum in HRS (hamc)", &pmom, 200, 1.01, 1.07);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "momb2", 
	   "Momentum in HRS (hamc)", &pmom, 200, 3.0, 3.17);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "mombI", 
	   "Momentum in HRS (hamc)", &pmom, 200, 2890, 3.180);


// Note, th0,ph0, etc are the initial values right after scattering.
// Keep in mind theta,phi here are actually tangents (transport notation)

   expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "thtrans",
		 "Init Theta transport",&th0,nbin,-0.3,0.3);
   expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "phtrans",
		 "Init Phi transport",&ph0,nbin,-0.3,0.3);
   expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "thph",
	  "Theta-Phi at target (unbiased)",&ph0,nbin,-0.3,0.3,
			  &th0,nbin,-0.3,0.3);
   expt->inout->BookHisto(kFALSE, kTRUE, ICOLLIM, "thphc",
	  "Theta-Phi at target (collimated)",&ph0,nbin,-0.3,0.3,
			  &th0,nbin,-0.3,0.3);
   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "thpha",
	  "Theta-Phi at target (accepted)",&ph0,nbin,-0.3,0.3,
			  &th0,nbin,-0.3,0.3);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "thphaw",
	  "Theta-Phi at target (accepted, weighted)",
                          &ph0,nbin,-0.3,0.3,
			  &th0,nbin,-0.3,0.3);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "thphawz",
	  "Theta-Phi at target (accepted, weighted, zoom)",
                          &ph0,nbin,-0.04,0.04,
			  &th0,nbin,-0.08,0.08);
   expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "dpp",
			  "dp/p generated ",
                          &dpp0,nbin,-0.007,0.002);

   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfpd",
			  "X-Y at focal plane detector",
			  &xfpd, nbin, -1, 1,
                          &yfpd, nbin, -0.2, 0.2);

// For designing the collimator
   expt->inout->BookHisto(kFALSE, kFALSE, ICOLLIM, "xycoll",
			  "X-Y at collimator",
			  &ytrans, nbin, -0.5, 0.5,
                          &xtrans, nbin, -0.5, 0.5);

   expt->inout->BookHisto(kFALSE, kTRUE, ICOLLIM, "xycolla",
			  "X-Y at collimator in acceptance",
			  &ytrans, nbin, -0.5, 0.5,
                          &xtrans, nbin, -0.5, 0.5);

// For designing the collimator2
   expt->inout->BookHisto(kFALSE, kFALSE, ICOLLIM2, "xycoll2",
			  "X-Y at collimator2",
			  &ytrans, nbin, -0.5, 0.5,
                          &xtrans, nbin, -0.5, 0.5);

   expt->inout->BookHisto(kFALSE, kTRUE, ICOLLIM2, "xycoll2a",
			  "X-Y at collimator2 in acceptance",
			  &ytrans, nbin, -0.5, 0.5,
                          &xtrans, nbin, -0.5, 0.5);



   htp = new TH2F("htp","Generated Theta-Phi",
              100,-0.3,0.3,100,-0.3,0.3);


   Float_t xbox = 0.8;
   Float_t ybox = 0.8;

// Note, since Xtrans is vertical, it makes a little more sense to 
// plot it on the vertical axis (and Ytrans on horizontal)

   expt->inout->BookHisto(kFALSE, kFALSE, ICOLLIM, "xycol", 
		      "Transport X-Y at collimator", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, ICOLLIM, "xycola", 
		   "Transport X-Y inside collimator acceptance", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);

   // But inside septum, X,Y are rotated 90, so put them back ...
   expt->inout->BookHisto(kFALSE, kFALSE, ISEPTIN, "xysepi", 
		      "Transport X-Y at Septum entrance", 
                            &xtrans, nbin,-xbox,xbox,
                            &ytrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, ISEPTIN, "xysepia", 
		   "Transport X-Y inside Sept-in acceptance", 
                            &xtrans, nbin,-xbox,xbox,
                            &ytrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, ISEPTIN, "xysepiaz", 
		   "Transport X-Y inside Sept-in accept (zoom)", 
                            &xtrans, nbin,-0.05,0.45,
                            &ytrans, nbin,-0.2,0.2);

   expt->inout->BookHisto(kFALSE, kFALSE, ISEPTOUT, "xysepo", 
		      "Transport X-Y at Septum exit", 
                            &xtrans, nbin,-xbox,xbox,
                            &ytrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, ISEPTOUT, "xysepoa", 
		   "Transport X-Y inside Sept-exit acceptance", 
                            &xtrans, nbin,-xbox,xbox,
                            &ytrans, nbin,-ybox,ybox);

   expt->inout->BookHisto(kFALSE, kFALSE, IQ1EXIT, "xyq1", 
		      "Transport X-Y at Q1 exit", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, IQ1EXIT, "xyq1a", 
		   "Transport X-Y inside Q1 acceptance", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);

   expt->inout->BookHisto(kFALSE, kFALSE, IDIPIN, "xydipi", 
		      "Transport X-Y at dipole entrance", 
 			    &ytrans, nbin,-1,1,
                            &xtrans, nbin,-6,-4);
   expt->inout->BookHisto(kFALSE, kTRUE, IDIPIN, "xydipia", 
		   "Transport X-Y inside dipole-in acceptance", 
			  &ytrans, nbin,-1,1,
			  &xtrans, nbin,-6,-4);

   expt->inout->BookHisto(kFALSE, kFALSE, IDIPEXIT, "xydipo", 
		      "Transport X-Y at dipole exit", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, IDIPEXIT, "xydipoa", 
		   "Transport X-Y inside dipole-exit acceptance", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);

   expt->inout->BookHisto(kFALSE, kFALSE, IQ3IN, "xyq3i", 
		      "Transport X-Y at Q3 entrance", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, IQ3IN, "xyq3ia", 
		   "Transport X-Y inside Q3-in acceptance", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);

   expt->inout->BookHisto(kFALSE, kFALSE, IQ3EXIT, "xyq3o", 
		      "Transport X-Y at Q3 exit", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);
   expt->inout->BookHisto(kFALSE, kTRUE, IQ3EXIT, "xyq3oa", 
		   "Transport X-Y inside Q3-exit acceptance", 
                            &ytrans, nbin,-xbox,xbox,
                            &xtrans, nbin,-ybox,ybox);

   expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "xyfoc1", 
		      "X-Y at focal plane (even if not accepted)", 
                            &ytrans, nbin,-0.5,0.5,
                            &xtrans, nbin,-1,1); 
   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "xyfoc2", 
		      "Accepted X-Y at focal plane", 
                            &ytrans, nbin,-0.1,0.1,
                            &xtrans, nbin,-0.6,0.2); 
    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfoc3", 
		      "Weighted X-Y at focal plane", 
                            &ytrans, nbin,-0.1,0.1,
                            &xtrans, nbin,-0.6,0.2); 

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfoc4", 
		      "Weighted X-Y at focal plane", 
                            &ytrans, nbin,-0.1,0.1,
                            &xtrans, nbin,-0.6,0.2); 

    expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "xyfoc5", 
		      "Unweighted X-Y at focal plane (X on X-axis)", 
                            &xtrans, nbin,-0.8,0.3,
                            &ytrans, nbin,-0.1,0.1);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfoc6", 
		      "Weighted X-Y at focal plane (X on X-axis)", 
                            &xtrans, nbin,-1,1,
                            &ytrans, nbin,-0.2,0.2);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfoc6z", 
		      "Weighted X-Y at focal plane (X on X-axis)", 
                            &xtrans, nbin,-0.1,0.08,
                            &ytrans, nbin,-0.08,0.07);


    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xfoc",
                         "X at focal plane",
			   &xtrans, nbin, -0.5, 0.1);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "yfoc",
                         "Y at focal plane",
			   &ytrans, nbin, -0.1, 0.1);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xdet",
                         "X (det. frame)",
			   &xdet, nbin, -1.2, 1.);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "qsq",
			   "Qsq (weighted, in accept)",
			   &qsq, 200,  0, 1);


    // Add some variables to the event ntuple
    // These get filled at the focal plane
    expt->inout->AddToNtuple("pmom",&pmom);
    expt->inout->AddToNtuple("dpp",&dpp);
    expt->inout->AddToNtuple("qsq",&qsq);
    expt->inout->AddToNtuple("theta",&theta);
    expt->inout->AddToNtuple("phi",&phi);
    expt->inout->AddToNtuple("xfoc",&xtrans);
    expt->inout->AddToNtuple("yfoc",&ytrans);
    expt->inout->AddToNtuple("xdet",&xdet);
    expt->inout->AddToNtuple("ydet",&ydet);
    expt->inout->AddToNtuple("mscol",&ms_collim);
    expt->inout->AddToNtuple("th0",&th0);
    expt->inout->AddToNtuple("xfpd",&xfpd);
    expt->inout->AddToNtuple("yfpd",&yfpd);
    expt->inout->AddToNtuple("thfpd",&thfpd);
    expt->inout->AddToNtuple("phfpd",&phfpd);
 

    did_init = kTRUE;

    return OK;      
}


Int_t hamcTrackOut::Generate(hamcExpt *expt) {

// To generate a track for an event.

  if ( did_init == kFALSE ) {
    cout << "hamcTrackOut::ERROR: uninitialized !"<<endl;
    return ERROR;
  }

  ms_collim = 0;  // reset
  thfpd = 0;
  phfpd = 0;
  xfpd  = 0;
  yfpd  = 0;

  hamcBeam *beam = expt->event->beam;

// Initial tvect X,Y,Z is tied to the beam
// X is vertical up, Y is horizontal.
// tvect units are meters.  
// (And the "angles" are actually tan(theta), tan(phi),
//    so those are approximately radians).

  if (beam) {
    Float_t xb,xt,yb,yt,zb,zt;
    xb = beam->GetTransX();  
    yb = beam->GetTransY();
    zb = beam->GetTransZ();
    xt = xb;
    yt = yb * TMath::Cos(theta_central);
    zt = zb + yb * TMath::Sin(theta_central);
    tvect->PutX(xt);
    tvect->PutY(yt);
    tvect->PutZ(zt);
  } else {
    tvect->Clear();
  }

  expt->physics->kine->Generate(expt);

  theta = expt->physics->kine->theta;
  phi = expt->physics->kine->phi;
  energy = expt->physics->kine->eprime;
  qsq = expt->physics->kine->qsq;

  if (energy < mass) energy = mass;  // extrema of rad tail

  pmom = TMath::Sqrt(energy*energy - mass*mass);
  
  dpp = 0;
  if (P0 != 0) dpp = (pmom - P0)/P0;

  tvect->PutDpp(dpp);

  LabToTrans();     // Get TRANSPORT angles.
  UpdateTrans();
  ComputePvect();   

  return OK;

}

Int_t hamcTrackOut::UpdateAtDet() {

  thfpd = tvect->GetTheta();
  phfpd = tvect->GetPhi();
  xfpd  = tvect->GetX() + det_dist*thfpd;
  yfpd  = tvect->GetY() + det_dist*phfpd;
  return OK;

}


Int_t hamcTrackOut::LabToTrans() {

// Convert lab angles to TRANSPORT.
// In lab, beam is along Z.
// theta, phi are the already-generated scattering angle
// and azimuthal angle in the lab frame.
// theta_central is the HRS optic axis angle with respect to Z.
// We need the transport vector (t) angles, actually need
// tan(theta_t), tan(phi_t) where 
// positive theta_t is vertical down
// phi_t sign depends on which spectrometer.
// for R-HRS positive phi_t is towards beam, and 
// for L-HRS positive phi_t is away from beam
 
  Float_t stc, ctc, ttc, sts, cts, sps, cps;
  stc = TMath::Sin(theta_central);
  ctc = TMath::Cos(theta_central);
  ttc = TMath::Tan(theta_central);
  sts = TMath::Sin(theta);
  cts = TMath::Cos(theta);
  sps = TMath::Sin(phi);
  cps = TMath::Cos(phi);

  Float_t uparam = stc*sps*sts + ctc*cts;

  if (uparam == 0) {
    cout << "hamcTrackOut:: ERROR: uparam = 0 !"<<endl;
    tvect->Clear();
    return ERROR;
  }

// The -1.0 is to obey sign convention, see above.
  Float_t tantheta_t = -1.0*cps*sts / uparam;

  Float_t tp_prime = sps*sts / uparam;

// Need to subtract the central scattering angle:
  Float_t xsign = 1.0;
  if (which_hrs == RIGHTHRS) xsign = -1.0;  // sign convention, see above
  Float_t tanphi_t = xsign*(tp_prime - ttc) / (1 + tp_prime*ttc);

  tvect->PutTheta(tantheta_t);
  tvect->PutPhi(tanphi_t);      

  th0 = tvect->GetTheta();
  ph0 = tvect->GetPhi();
  x0  = tvect->GetX();
  y0  = tvect->GetY();
  z0  = tvect->GetZ();
  dpp0 = tvect->GetDpp();

  htp->Fill(tanphi_t, tantheta_t);

  *tvect_orig = *tvect;  // update the "origin" tvect.

  //  cout << "stuff "<<uparam<<"  "<<tantheta_t<<"  "<<tp_prime<<"  "<<tanphi_t<<endl;

  if (debug) {
    cout << "LabToTrans print"<<endl;
    tvect->Print();
  }

  return OK;

}

void hamcTrackOut::ComputePvect() {
  
// Here, Z is along beam. 
// X is to the right in horizontal plane.
// Y is vertical.

  Float_t sts, cts, sps, cps;
  sts = TMath::Sin(theta);
  cts = TMath::Cos(theta);
  sps = TMath::Sin(phi);
  cps = TMath::Cos(phi);

  plab_x = sps*sts * pmom;
  plab_y = cps*sts * pmom;
  plab_z = cts * pmom;

}

void hamcTrackOut::ComputeQsqToTrack(const hamcTrack *trk) {

// Could compute Qsq of *this to *trk, but dont really need
// this method because hamcKine does it for you.

  qsq = -9999;        
  if (!trk) return;

  Float_t e2,px2,py2,pz2;

  e2  = trk->GetEnergy();
  px2 = trk->GetPx();
  py2 = trk->GetPy();
  pz2 = trk->GetPz();

  qsq = (energy - e2)*(energy-e2) - (
        (plab_x - px2)*(plab_x - px2) + 
        (plab_y - py2)*(plab_y - py2) + 
        (plab_z - pz2)*(plab_z - pz2) );

  qsq = -1*qsq;  // make it positive

}



   
