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
#include "hamcTransGuido.h"
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
  det_dist= 0; // meters
  trktype = "out";
  tvect_guido = 0;
  xgui = 0; ygui = 0; thgui = 0; phgui = 0; 
}


hamcTrackOut::hamcTrackOut(string pid, Float_t ee, Float_t x, Float_t th, Float_t y, Float_t ph, Float_t dp) : hamcTrack(pid,ee,x,th,y,ph,dp)
{
  did_init = kFALSE;
  trktype = "out";
  which_hrs = 0;
  tvect_guido = 0;
  xgui = 0; ygui = 0; thgui = 0; phgui = 0; 
}

hamcTrackOut::~hamcTrackOut() {
  if (tvect_guido) delete tvect_guido;
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
  parser.Load(expt->inout->GetStrVect("hrs_setup"));
  Int_t use_guido=0;
  if (parser.IsFound("useguido")) {
    // Would rather do 'if (spec->IsGuidoTrans())' but 
    // spec is not initialized yet.
     use_guido=1;
     tvect_guido = new hamcTransVect();
  }

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
	        "Theta Generated", &theta, nbin,0.8*thetamin,1.2*thetamax);
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
                &theta, 200, 0.05,0.15);
   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "th5", 
		"Theta at focal plane (unweighted, in accept)", 
                &theta, nbin, 0.8*thetamin,1.2*thetamax);

// To measure acceptance function
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "hthin", 
	        "Theta Generated", &theta_deg, 600,2.0,8.0);
   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "hthacc", 
	        "Theta Accepted", &theta_deg, 600,2.0,8.0);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "hthaccw", 
	        "Theta Accepted, Weighted", &theta_deg, 600,2.0,8.0);
   expt->inout->BookHisto(kTRUE, kTRUE, ITARGET, "xtg","X at target",&xtrans,200,-0.1,0.1);
   expt->inout->BookHisto(kTRUE, kTRUE, ITARGET, "ytg","Y at target",&ytrans,200,-0.05,0.05);
   expt->inout->BookHisto(kFALSE, kFALSE, ITARGET, "ph", 
		"Phi at target", &phi, nbin, -0.2,1.2*phimax);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "mom", 
	        "Momentum in HRS", &pmom, nbin, 0.5*P0,1.05*P0);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "momu", 
	        "Momentum in HRS", &pmom, nbin, 0.5*P0,1.1*P0);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "momz", 
		 "Momentum in HRS (zoom)", &pmom, 100, 0.925, 0.955);
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
	  "Theta-Phi at ztarget (accepted)",&ph0,nbin,-0.032,0.032,
			  &th0,nbin,-0.067,0.067);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "thphaw",
	  "Theta-Phi at target (accepted, weighted)",
                          &ph0,nbin,-0.1,0.1,
			  &th0,nbin,-0.1,0.1);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "thphawz",
	  "Theta-Phi at target (accepted, weighted, zoom)",
                          &ph0,100,-0.032,0.032,
			  &th0,100,-0.067,0.067);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "thphawzms",
	  "Theta-Phi at target (accepted, weighted, zoom, mult scatt)",
                          &phtgt,200,-0.032,0.032,
			  &thtgt,200,-0.067,0.067);
//                          &ph_ms,200,-0.032,0.032,
//			  &th_ms,200,-0.067,0.067);

   expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "dpp",
			  "dp/p generated ",
                          &dpp0,nbin,-0.007,0.002);

   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfpd",
			  "X-Y at focal plane detector",
			  &xfpd, nbin, -1, 1,
                          &yfpd, nbin, -0.2, 0.2);

   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "th0",
                          "Theta(transport) at target",
                          &thtgt, 100, -0.06, 0.06);

   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "phi0",
                          "Phi(transport) at target",
                          &ph0, 100, -0.04, 0.04);

   expt->inout->BookHisto(kTRUE, kTRUE, ICOLLIM, "phtrans1",
                          "Phi(transport) with mult. scatt",
                          &phtrans, nbin, -0.04, 0.04);

   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "phtrans2",
                          "Phi(transport) with mult. scatt",
                          &phtrans, nbin, -0.04, 0.04);


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

   // Plots to compare to data
  expt->inout->BookHisto(kFALSE, kFALSE, IFOCAL, "thph2",
		     "Theta-Phi at target (Monte Carlo)",&ph0,nbin,-0.06,0.06,
			  &th0,nbin,-0.1,0.1);
   expt->inout->BookHisto(kFALSE, kTRUE, ICOLLIM, "thphc2",
		  "Theta-Phi at target (collimated)",&ph0,nbin,-0.06,0.06,
			  &th0,nbin,-0.1,0.1);
   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "thpha2",
	  "Theta-Phi at target (accepted)",&ph0,nbin,-0.06,0.06,
			  &th0,nbin,-0.1,0.1);
   expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "thpha4",
	  "Theta-Phi at target (accepted, weighted)",&ph0,nbin,-0.06,0.06,
			  &th0,nbin,-0.1,0.1);

// For designing the collimator3
   expt->inout->BookHisto(kFALSE, kFALSE, ICOLLIM3, "xycoll3",
			  "X-Y at collimator3",
			  &ytrans, nbin, -0.5, 0.5,
                          &xtrans, nbin, -0.5, 0.5);

   expt->inout->BookHisto(kFALSE, kTRUE, ICOLLIM3, "xycoll3a",
			  "X-Y at collimator3 in acceptance",
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
			  &ytrans, nbin,-0.2,0.2,
			  &xtrans, nbin,-0.35,-0.17);

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
                            &ytrans, nbin,-0.8,0.8,
                            &xtrans, nbin,-0.8,0.8);

   expt->inout->BookHisto(kFALSE, kTRUE, IQ1EXIT, "xyq1a", 
		   "Transport X-Y inside Q1 acceptance", 
                            &ytrans, nbin,-0.8,0.8,
                            &xtrans, nbin,-0.8,0.8);

   Float_t xdlo,xdhi;
   parser.Load(expt->inout->GetStrVect("hrs_setup"));
   if (parser.IsFound("thrstrans") || parser.IsFound("hrstrans")) {
     cout << "hamcTrackOut:: Is THRSTrans "<<endl;
     xdlo = -1; xdhi = 1;
   } else {
     cout << "hamcTrackOut:: Is NOT THRSTrans "<<endl;
     xdlo = -5.5; xdhi = -4.8;
   }

   expt->inout->BookHisto(kFALSE, kFALSE, IDIPIN, "xydipi", 
		      "Transport X-Y at dipole entrance", 
			  &ytrans, nbin,-0.4,0.4,
			  &xtrans, nbin,xdlo,xdhi);

   expt->inout->BookHisto(kFALSE, kTRUE, IDIPIN, "xydipia", 
		   "Transport X-Y inside dipole-in acceptance", 
			  &ytrans, nbin,-0.4,0.4,
			  &xtrans, nbin,xdlo,xdhi);

   expt->inout->BookHisto(kFALSE, kFALSE, IDIPIN, "xydipsi", 
		      "(std) Trans X-Y at dipole entrance", 
 			    &ytrans, nbin,-1,1,
                            &xtrans, nbin,-1,1);

   expt->inout->BookHisto(kFALSE, kTRUE, IDIPIN, "xydipsia", 
		   "(std) Trans X-Y inside dipole-in acceptance", 
			  &ytrans, nbin,-1,1,
			  &xtrans, nbin,-1,1);

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
			  &xtrans, nbin,-1.2,0.2,
			  &ytrans, nbin,-0.4,0.4);

   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "xyfoc1a", 
		      "Accepted X-Y at focal plane", 
			  &xtrans, nbin,-1.2,0.2,
			  &ytrans, nbin,-0.4,0.4);

   expt->inout->BookHisto(kFALSE, kTRUE, IFOCAL, "xyfoc2a", 
		      "Accepted X-Y at focal plane", 
			  &xtrans, 100,-0.8,0.3,
			  &ytrans, 100,-0.05,0.05);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfoc2", 
		      "Weighted  X-Y at focal plane", 
                            &xtrans, 200,-1.2,0.6,
 			    &ytrans, 200,-0.2,0.2);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xyfoc3", 
		      "Weighted  X-Y at focal plane", 
                            &xtrans, 200,-1,0.25,
                            &ytrans, 200,-0.1,0.1); 

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xfoc",
                         "X at focal plane",
			   &xtrans, 200, -1.2, 0.6);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "yfoc",
                         "Y at focal plane",
			   &ytrans, 200, -0.05, 0.05);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "xdet",
                         "X (det. frame)",
			   &xdet, nbin, -1.2, 1.);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "qsq",
			   "Qsq (weighted, in accept)",
			   &qsq, 100,  0.003, 0.014);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "qsq_hap3",
			   "Qsq (weighted, in accept)",
			   &qsq, 200,  0.3, 0.8);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "qsq_obs",
			   "Qsq (weighted, in accept)",
			   &qsq_obs, 200, 0.0015, 0.025);

    expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "qsq_atrk",
			   "Qsq (weighted, in accept)",
			   &qsq_atrk, 200, 0.0015, 0.025);


    // Add some variables to the event ntuple
    // These get filled at the focal plane

    expt->inout->AddToNtuple("pmom",&pmom);
    expt->inout->AddToNtuple("dpp",&dpp);
    expt->inout->AddToNtuple("qsq",&qsq);
    expt->inout->AddToNtuple("qsq_obs",&qsq_obs);
    expt->inout->AddToNtuple("qsq_atrk",&qsq_atrk);
    expt->inout->AddToNtuple("qsqfr",&qsqfr);
    expt->inout->AddToNtuple("theta",&theta);
    expt->inout->AddToNtuple("phi",&phi);
    expt->inout->AddToNtuple("xfoc",&xtrans);
    expt->inout->AddToNtuple("yfoc",&ytrans);
    expt->inout->AddToNtuple("mscol",&ms_collim);
    expt->inout->AddToNtuple("th0",&th0);
    expt->inout->AddToNtuple("ph0",&ph0);
    expt->inout->AddToNtuple("xfpd",&xfpd);
    expt->inout->AddToNtuple("yfpd",&yfpd);
    expt->inout->AddToNtuple("thfpd",&thfpd);
    expt->inout->AddToNtuple("phfpd",&phfpd);
// at tgt, after mult. scatt. and resol. smearing -- to compare to data
    expt->inout->AddToNtuple("thtgt",&thtgt);  
    expt->inout->AddToNtuple("phtgt",&phtgt);
    
    expt->inout->AddToNtuple("thchk",&thchk);
    expt->inout->AddToNtuple("phms",&ph_ms);
    expt->inout->AddToNtuple("thms",&th_ms);

    if (use_guido) {
      expt->inout->AddToNtuple("xgui",&xgui);
      expt->inout->AddToNtuple("ygui",&ygui);
      expt->inout->AddToNtuple("thgui",&thgui);
      expt->inout->AddToNtuple("phgui",&phgui);
    }

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

  Int_t status =  expt->physics->kine->Generate(expt);

  if (status == -1) {
    return status;
  }

  theta = expt->physics->kine->theta;
  phi = expt->physics->kine->phi;
  theta_deg = theta * 180.0/PI;
  phi_deg = phi * 180.0/PI;
  energy = expt->physics->kine->eprime;
  qsq = expt->physics->kine->qsq;

  if (energy < mass) energy = mass;  // extrema of rad tail

  pmom = TMath::Sqrt(energy*energy - mass*mass);
  
  dpp = 0;
  if (P0 != 0) dpp = (pmom - P0)/P0;

  tvect->PutDpp(dpp);

  LabToTrans(expt);     // Get TRANSPORT position and angles.
  UpdateTrans();
  ComputePvect();   

  return OK;

}

Int_t hamcTrackOut::GenerateOut(hamcExpt *expt) {

// It's assumed that multiple scattering in the target
// has now been applied in hamcEvent

  if ( did_init == kFALSE ) {
    cout << "hamcTrackOut::ERROR: uninitialized !"<<endl;
    return ERROR;
  }

// Generate the outgoing track (after scattering), i
// including "observed" Qsq = qsq_obs.  

  expt->physics->kine->GenerateOut(expt);

  qsq_obs = expt->physics->kine->qsq_obs;
  qsq_atrk = expt->physics->kine->qsq_atrk;
  qsqfr = expt->physics->kine->qsqfr;

  return OK;

}

Int_t hamcTrackOut::UpdateAtDet() {

  thfpd = tvect->GetTheta();
  phfpd = tvect->GetPhi();
  xfpd  = tvect->GetX() + det_dist*thfpd;
  yfpd  = tvect->GetY() + det_dist*phfpd;
  return OK;

}

Int_t hamcTrackOut::UpdateGuidoFocal(hamcSpecHRS *spec) {

  // Initialize Guido's transport vector
  if (!tvect_guido) {
    cout << "hamcTrackOut::Error: No Guido tvect "<<endl;
    return 0;
  }

  *tvect_guido = *tvect_orig;

  // Transform the vector 
  if (spec->IsGuidoTrans()) {
    spec->tguido->TransForm(tvect_guido,IFOCAL);
  }
     
  xgui  = tvect_guido->GetX();
  ygui  = tvect_guido->GetY();
  thgui = tvect_guido->GetTheta();
  phgui = tvect_guido->GetPhi();

  return OK;

}

Int_t hamcTrackOut::LabToTrans(hamcExpt *expt) {

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

  Int_t ldebug=0;
  Int_t useold_test = 0;

  Float_t stc, ctc, ttc, sts, cts, tts, sps, cps;
  stc = TMath::Sin(theta_central);
  ctc = TMath::Cos(theta_central);
  ttc = TMath::Tan(theta_central);
  sts = TMath::Sin(theta);
  cts = TMath::Cos(theta);
  tts = TMath::Tan(theta);
  sps = TMath::Sin(phi);
  cps = TMath::Cos(phi);

  hamcBeam *beam = expt->event->beam;

// Initial tvect X,Y,Z is tied to the beam
// X is vertical up, Y is horizontal, positive to right.
// tvect units are meters.  
// (And the "angles" are actually tan(theta), tan(phi),
//    so those are approximately radians).

  xsign = 1.0;
  if (which_hrs == LEFTHRS) xsign = -1.0;  // sign convention, see above

  if (beam) {
    Float_t xb,xt,yb,yt,zb,zt;
    xb = beam->GetTransX();  
    yb = beam->GetTransY();
    zb = beam->GetTransZ();
    xt = xb  - (yb * sts * cps * stc + zb * sts * cps * ctc)/(cts * ctc + sts * sps * stc);
    yt = ( yb - zb * tts * sps * sps)/(ctc + tts * sps * stc);
    zt = zb;  // actually irrelevant

    if (useold_test) {
       xt = xb;
       yt = yb * TMath::Cos(theta_central);
       zt = zb * TMath::Cos(theta_central) + yb * TMath::Sin(theta_central); 
    }

    tvect->PutX(xt);
    tvect->PutY(yt);
    tvect->PutZ(zt);
  } else {
    tvect->Clear();
  }

  Int_t flip_phi = 0;
  if (expt->GetSpectrom(0)->IsWarmSeptum()) {
// Use R-HRS for transport, but flip ph0 for comparing to real data
    if (which_hrs == LEFTHRS) flip_phi = 1;
    xsign = 1.0;  
  }

// The following 3 lines are from D. Wang (new, June 2009)
// tanphprime is slightly different and correct now.

  Float_t tanphprime = sts*sps/cts;
  Float_t tanphi_t = xsign * ( ttc - tanphprime)/(1 + ttc * tanphprime); 
  Float_t tantheta_t = -1. *  (sts * cps) * TMath::Sqrt( 1 + tanphi_t * tanphi_t)/(TMath::Sqrt(cts * cts + sts*sts*sps*sps));

  ph0 = tanphi_t;
  if (flip_phi) ph0 = -1.0*tanphi_t;
  th0 = tantheta_t;
  ph_ms = tanphi_t;
  th_ms = tantheta_t;

  tvect->PutTheta(tantheta_t);
  tvect->PutPhi(tanphi_t);      

// Check of scattering angle
  thchk = TMath::ACos(( ctc + xsign*fabs(stc)*tanphi_t)/(TMath::Sqrt(1 + tanphi_t * tanphi_t + tantheta_t * tantheta_t)));

  //  cout << "Origin theta, phi "<<th0<<"   "<<ph0<<endl;
 
  x0  = tvect->GetX();
  y0  = tvect->GetY();
  z0  = tvect->GetZ();
  dpp0 = tvect->GetDpp();

  htp->Fill(tanphi_t, tantheta_t);

  Float_t th0resol = 0.0;
  Float_t ph0resol = 0.0;

  th0sm = th0 + th0resol*gRandom->Gaus();
  ph0sm = ph0 + ph0resol*gRandom->Gaus();

  *tvect_orig = *tvect;  // update the "origin" tvect.

  if (ldebug) {
    cout << "\n\nDebug LabToTrans for whichhrs ";
    cout <<  which_hrs<<"  "<<RIGHTHRS<<endl;
    cout << "Transport angles "<<endl;
    cout << "tan(theta) = "<<tantheta_t<<endl;
    cout << "tan(phi)   = "<<tanphi_t<<endl;
    cout << "Scat angle chk = "<<theta<<"  "<<thchk<<endl;
    cout << "Scat angle (degr) = "<<180*theta/PI<<"  ";
    cout << 180*thchk/PI<<endl;
    cout << endl<<"LabToTrans print"<<endl;
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

  pnoms_x = sps*sts * pmom;
  pnoms_y = cps*sts * pmom;
  pnoms_z = cts * pmom;

// Initially, plab is the momentum without multiple scattering.
// This may be update later after Mult. Scatt. in
// hamcTrack::MultScatt called by the event loop in hamcEvent 

  plab_x = pnoms_x;   
  plab_y = pnoms_y;
  plab_z = pnoms_z;

}

void hamcTrackOut::ComputeQsqToTrack(const hamcTrack *trk) {

// Obsolete method.
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


Float_t hamcTrackOut::GetUpdatedScatt() {

// Get, on demand, the updated scattering angle constructed from tg_th and tg_ph
// This will include whatever multiple scattering and resolution effects have occured.

  Float_t stc, ctc, tantheta_t, tanphi_t, thloc;

  stc = TMath::Sin(theta_central);
  ctc = TMath::Cos(theta_central);

  tantheta_t = tvect->GetTheta();
  tanphi_t = tvect->GetPhi();

// Reconstructed scattering angle

  thloc = TMath::ACos(( ctc + xsign*fabs(stc)*tanphi_t)/(TMath::Sqrt(1 + tanphi_t * tanphi_t + tantheta_t * tantheta_t)));

  return thloc;

}
   


   
