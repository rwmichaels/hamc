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

  expt->inout->AddToNtuple("inaccept",&inaccept);
  expt->inout->AddToNtuple("xcol",&xcol);
  expt->inout->AddToNtuple("ycol",&ycol);
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

// Weight initially 1 (used by histograms)
   expt->inout->SetWeight(1);  

   expt->target->Zscatt();

   if (expt->physics->Radiate(expt) == -1) return 1;

   if (beam)
     beam->Generate(expt);

   inaccept = 1;
   xcol = -999;  ycol = -999;
   xsep = -999;  ysep = -999;

// Loop over spectrometers

   for (Int_t ispect=0; ispect<expt->GetNumSpectrom(); ispect++) {

     hamcSpecHRS *spect = expt->GetSpectrom(ispect);
     hamcTrackOut *track = trackout[ispect]; // Assumes 1 track per spectrom.
     if (!track) continue;

     if (track->Generate(expt)==-1) {
       //       cout<<"no good dis event generated"<<endl;
       return OK;
     }

     expt->physics->Generate(expt); 


// Weight by cross section (optionally used for some histograms)
     Float_t weight = 1e5*expt->physics->GetCrossSection();
 
     expt->inout->SetWeight(weight);

     track->MultScatt(expt, ITARGET);

     Float_t xb4trans = track->GetTransX();
     Float_t yb4trans = track->GetTransY();
     //    cout<<"xb4trans="<<xb4trans<<"yb4trans="<<yb4trans<<endl;
     //     track->xyb4trans->Fill(yb4trans, xb4trans);

     Float_t tb4trans = track->GetTransTheta();
     Float_t pb4trans = track->GetTransPhi();
     //     track->tpb4trans->Fill(tb4trans,pb4trans);

     Float_t Dpb4trans = track->GetTransDp();
     //     track->dpb4trans->Fill(Dpb4trans);


// Loop over break points in spectrometer

     if (debug) {
      cout << "==================================="<<endl;
      cout << "ONE EVENT, starting break points "<<endl;
     }

     for (Int_t ibrk=0; ibrk<spect->GetNumBrk(); ibrk++) {

       brkpoint = spect->break_point[ibrk]->where;

       if (brkpoint != ITARGET) {

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
         if (brkpoint == ICOLLIM || brkpoint == ICOLLIM2) {
	   xcol = track->GetTransX();
	   ycol = track->GetTransY();
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

