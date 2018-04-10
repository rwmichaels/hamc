//  hamcTransLer4deg   -- Transport model using LeRose transfer functions
//  for the 4 degree septum for CREX 
//  R. Michaels  April 2013

// How many underscores ?  
// use "nm" to tell, "nm crex_4degr.o" and if check if you see one or two 
// trailing underscores


#include "hamcTransLer4deg.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTransLer4deg)
#endif


hamcTransLer4deg::hamcTransLer4deg() 
{
}

hamcTransLer4deg::~hamcTransLer4deg() {
}

Int_t hamcTransLer4deg::Init(hamcSpecHRS *spect) {

  which_spectrom = spect->which_spectrom;
  collim_distance = spect->GetCollimDist();

  return OK;

}

void hamcTransLer4deg::Print() {
}

Int_t hamcTransLer4deg::TransForm(hamcTrack *trk, Int_t where) const {

// Transforms a track "trk" to the destination = where.
// The origin of the track is known (trk->origin = ITARGET, ICOLLIM, etc).
// The LeRose functions are used.
// They are only defined for certain combinations of (origin,where).
// Any undefined combination of (origin,where) will be ignored,
// i.e. no change to the track.

// Input is transport vector. +X is down.
// The output is in a 'funny' Snake coordinate system.
// Therefore you should not use the output as a launch for a new track
// unless you make the proper coordinate transformation; however, 
// the (x,y) can be compared to the aperture function which does
// indeed take into account the new coordinates (this is why, for
// example, some apertures appear not centered in the coordinates).
// The final output at the focal plane is, once again, transport.

// NOTE: For warm septum there is no distinction between Left
// and Right HRS at the moment.


  Int_t origin = trk->origin;
  Int_t xbig = -9999;
  Int_t debug_trans = 1; // to get do some debugging here.

  int dimen=5;
  float xtrans[5], xout[5];

  xtrans[0] = trk->tvect_orig->GetX();
  xtrans[1] = trk->tvect_orig->GetTheta();
  xtrans[2] = trk->tvect_orig->GetY();
  xtrans[3] = trk->tvect_orig->GetPhi();
  xtrans[4] = trk->tvect_orig->GetDpp();

  for (Int_t i=0; i<4; i++) xout[i]=xbig;
  xout[4] = xtrans[4];

   if (origin == ITARGET) {

    switch(where) {

      case ISEPTIN: 
        break;

      case ISEPTOUT: 
 
        xout[0]=x_s4_sext_(xtrans,&dimen);
        xout[1]=t_s4_sext_(xtrans,&dimen);
        xout[2]=y_s4_sext_(xtrans,&dimen);
        xout[3]=p_s4_sext_(xtrans,&dimen);

        break;

      case ICOLLIM3: 

	trk->tvect->Load(xtrans);
        Drift(0.797, trk);  // 79.7 cm to sieve.
        return OK;

      case ICOLLIM: 
      case ICOLLIM2: 

        xout[0]=x_s4_q1en_(xtrans,&dimen);
        xout[1]=t_s4_q1en_(xtrans,&dimen);
        xout[2]=y_s4_q1en_(xtrans,&dimen);
        xout[3]=p_s4_q1en_(xtrans,&dimen);
        break;

      case IQ1EXIT: 

        xout[0]=x_s4_q1ex_(xtrans,&dimen);
        xout[1]=t_s4_q1ex_(xtrans,&dimen);
        xout[2]=y_s4_q1ex_(xtrans,&dimen);
        xout[3]=p_s4_q1ex_(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_s4_den_(xtrans,&dimen);
        xout[1]=t_s4_den_(xtrans,&dimen);
        xout[2]=y_s4_den_(xtrans,&dimen);
        xout[3]=p_s4_den_(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_s4_dex_(xtrans,&dimen);
        xout[1]=t_s4_dex_(xtrans,&dimen);
        xout[2]=y_s4_dex_(xtrans,&dimen);
        xout[3]=p_s4_dex_(xtrans,&dimen);

        break;

      case IQ3IN:

        xout[0]=x_s4_q3en_(xtrans,&dimen);
        xout[1]=t_s4_q3en_(xtrans,&dimen);
        xout[2]=y_s4_q3en_(xtrans,&dimen);
        xout[3]=p_s4_q3en_(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_s4_q3ex_(xtrans,&dimen);
        xout[1]=t_s4_q3ex_(xtrans,&dimen);
        xout[2]=y_s4_q3ex_(xtrans,&dimen);
        xout[3]=p_s4_q3ex_(xtrans,&dimen);
        break;

      case IFOCAL:

        xout[0]=x_s4_fp_(xtrans,&dimen);
        xout[1]=t_s4_fp_(xtrans,&dimen);
        xout[2]=y_s4_fp_(xtrans,&dimen);
        xout[3]=p_s4_fp_(xtrans,&dimen);
        break;

      case IPLANE1:

	xout[0]=x_s4_fp_(xtrans,&dimen);
        xout[1]=t_s4_fp_(xtrans,&dimen);
        xout[2]=y_s4_fp_(xtrans,&dimen);
        xout[3]=p_s4_fp_(xtrans,&dimen);

        trk->tvect->Load(xout);
        Drift(1.0, trk);  // 1 meter
        return OK;

      case IPLANE2:

	xout[0]=x_s4_fp_(xtrans,&dimen);
        xout[1]=t_s4_fp_(xtrans,&dimen);
        xout[2]=y_s4_fp_(xtrans,&dimen);
        xout[3]=p_s4_fp_(xtrans,&dimen);

        trk->tvect->Load(xout);
        Drift(1.48, trk);  // 1.48 m
        return OK;

      default:

        cout << "TransLer::WARNING: undefined destination "<<where<<endl;

    }

  }
  
  if (origin == ICOLLIM || origin == ICOLLIM2 || origin == ICOLLIM3) {

        cout << "TransLer4deg::WARNING: Cannot use collimator origin "<<origin<<endl;

  }

  if (xout[0] != xbig) { // this means xout was loaded; if it wasn't
                         // loaded then we dont change trk->tvect.

     trk->tvect->Load(xout);

  }

  
  return OK;

}

