//  HAMCEVENT   -- One events in hamc
//  R. Michaels  June 2008

#include "hamcEvent.h"
#include "hamcExpt.h"
#include "hamcSpecHRS.h"
#include "hamcPhysics.h"
#include "hamcBeam.h"
#include "hamcKine.h"
#include "hamcTrack.h"
#include "hamcTrackOut.h"
#include "hamcInout.h"
#include "hamcTarget.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcEvent)
#endif
 
hamcEvent::hamcEvent() 
{
  evnum = 0;
  beam = new hamcBeam();
  did_init = kFALSE;
}

hamcEvent::~hamcEvent()  
{
  if (beam) delete beam;
  for (std::vector<hamcTrackOut*>::iterator ith = trackout.begin();
      ith != trackout.end(); ith++) delete *ith;
}

Int_t hamcEvent::InitBeam(hamcExpt* expt) {

  if (beam) beam->Init(expt);
  return OK;

}

Int_t hamcEvent::Init(hamcExpt* expt) {

  Int_t nspect = expt->GetNumSpectrom();

  if (nspect == 0) {
    cout << "hamcEvent::ERROR: no spectrometers ?!"<<endl;
  } else {
    for (Int_t i=0; i<nspect; i++) {
       trackout.push_back(new hamcTrackOut());
       trackout[i]->Init(i,expt);
    }
  }

  tracknoms = new hamcTrackOut();
  tracksieve = new hamcTrackOut();
  tracksievenoms = new hamcTrackOut();

  expt->inout->AddToNtuple("inaccept",&inaccept);
  expt->inout->AddToNtuple("xtgt",&xtgt);
  expt->inout->AddToNtuple("ytgt",&ytgt);
  expt->inout->AddToNtuple("thtgt",&thtgt);
  expt->inout->AddToNtuple("phtgt",&phtgt);
  expt->inout->AddToNtuple("xcol",&xcol);
  expt->inout->AddToNtuple("ycol",&ycol);
  expt->inout->AddToNtuple("thcol",&thcol);
  expt->inout->AddToNtuple("phcol",&phcol);
  expt->inout->AddToNtuple("thnoms",&thnoms);
  expt->inout->AddToNtuple("phnoms",&phnoms);
  expt->inout->AddToNtuple("thsiv",&thsiv);
  expt->inout->AddToNtuple("phsiv",&phsiv);
  expt->inout->AddToNtuple("thsinoms",&thsinoms);
  expt->inout->AddToNtuple("phsinoms",&phsinoms);
  expt->inout->AddToNtuple("xsep",&xsep);
  expt->inout->AddToNtuple("ysep",&ysep);

  did_init = kTRUE;

  return OK;

}

Int_t hamcEvent::Process(hamcExpt* expt) {

// To process this one event 

   if ( did_init == kFALSE ) {
     cout << "hamcEvent::WARNING: did not initialize the event class"<<endl;
     return OK;
   }

   evnum++;

// Weight initially 1 (used by histograms)
   expt->inout->SetWeight(1);  

   expt->target->Zscatt();

   if (expt->physics->Radiate(expt) == -1) return 1;

   if (beam) beam->Generate(expt);

   inaccept = 1;
   xtgt = -999;  ytgt = -999;
   thtgt = -999; phtgt = -999;
   xcol = -999;  ycol = -999;
   thcol = -999; phcol = -999;
   thnoms = -999; phnoms = -999;
   thsiv = -999; phsiv = -999;
   thsinoms = -999; phsinoms = -999;
   xsep = -999;  ysep = -999;

// Loop over spectrometers

   for (Int_t ispect=0; ispect<expt->GetNumSpectrom(); ispect++) {

     hamcSpecHRS *spect = expt->GetSpectrom(ispect);
     hamcTrackOut *track = trackout[ispect]; // Assumes 1 track per spectrom.
     
     if (!track) continue;

     Int_t trkstat = track->Generate(expt);

     //     cout << "\n\nInitial track out "<<endl;
     //     track->Print();

     *tracknoms = *track;
     *tracksievenoms = *track;

     track->MultScatt(expt, ITARGET_FULL);

     *tracksieve = *track;

     expt->physics->Generate(expt); 

     track->GenerateOut(expt);

     // used in the breakpoint loop, after target
     Float_t angcut = expt->GetAngCut();

// Weight by cross section (optionally used for some histograms)
// The 100 is just for convenience.

     Float_t weight = 100.*expt->physics->GetCrossSection(); 
 
     expt->inout->SetWeight(weight);

     if (trkstat == -1) {  // but kill event if track not generated
       expt->inout->SetWeight(0);  
     }

// Loop over break points in spectrometer

     if (debug) {
      cout << "==================================="<<endl;
      cout << "ONE EVENT, starting break points "<<endl;
     }

     for (Int_t ibrk=0; ibrk<spect->GetNumBrk(); ibrk++) {

       brkpoint = spect->break_point[ibrk]->where;

       if (brkpoint == ITARGET) {
	   xtgt = track->GetTransX();
	   ytgt = track->GetTransY();
           thtgt = track->GetTransTheta();
           phtgt = track->GetTransPhi();
       }

       if (brkpoint > ITARGET && angcut > 0) {  // To simulate mistuned Septum need 
          Float_t theta_loc = (180./PI)*track->GetUpdatedScatt();
           if (theta_loc < angcut) {
             inaccept = 0;
           }
           if ((evnum%10000)==0) cout << "Warning !  Using (w/MS) angle cut "<<angcut<<endl;
        }



       if (brkpoint != ITARGET) {

         tracknoms->Transport(spect->transport, brkpoint);
         track->Transport(spect->transport, brkpoint);

// Here it's assumed the break points are exclusive
// (2 points might be at same Z but don't overlap, so the track
// is in one or the other).  Note, if we MultScatt (it only
// happens if radlen !=0) we reset the origin of the track.


         if ( track->InAccept(spect->Aperture(ibrk))  )  { // Track is in acceptance

           track->MultScatt(spect->Aperture(ibrk), brkpoint);

         }  else {  

           if (spect->break_point[ibrk]->IsCut()) {   // If we care about acceptance
                   inaccept = 0;
	   }
	 }

// Things to do for specific break points
         if (brkpoint == ICOLLIM || brkpoint == ICOLLIM2 || brkpoint == ICOLLIM3|| brkpoint == ICOLLIM4 ) {
	   xcol = track->GetTransX();
	   ycol = track->GetTransY();
           thcol = track->GetTransTheta();
           phcol = track->GetTransPhi();
           thnoms = tracknoms->GetTransTheta();
           phnoms = tracknoms->GetTransPhi();
           thsiv = tracksieve->GetTransTheta();
           phsiv = tracksieve->GetTransPhi();
           thsinoms = tracksievenoms->GetTransTheta();
           phsinoms = tracksievenoms->GetTransPhi();
           Float_t xcut = 0.043; 
           Float_t ycut = 0.0175;
           Float_t zzz = 0.83;
           Float_t xchk = zzz*track->th0;
	   Float_t ychk = zzz*track->ph0;
	   // This is just a check.  You may also do this in the ntuple.
	   //     cout << "Extrap. from tgt:  "<<zzz<<"  "<<track->th0<<"  "<<track->ph0<<"  "<<xchk<<"  "<<ychk<<endl;
	   //           if (xchk < -1*xcut || xchk > xcut) inaccept = 0;
	   //           if (ychk < -1*ycut || ychk > ycut) inaccept = 0;
 	 }
         if (brkpoint == ICOLLIM2) {
	     track->Eloss(expt, spect->Aperture(ibrk), brkpoint);
      	 }
         if (brkpoint == ISEPTIN) {
	   xsep = track->GetTransX();
	   ysep = track->GetTransY();
	 }
         if (brkpoint == IFOCAL) {
           if (spect->IsGuidoTrans()) {
	      track->UpdateGuidoFocal(spect);
	   }
           track->UpdateAtDet();
	 }
       }
	     
// The following line fills histograms, and will fill
// the ntuple at the focal plane when brkpoint=IFOCAL.

       expt->inout->Process(expt);

     }  // loop over break points

   }    // spectrometer loop

   expt->EventAnalysis();

   return OK;
  
}

