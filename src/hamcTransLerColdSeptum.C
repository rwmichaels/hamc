//  hamcTransLerColdSeptum   -- Transport model using LeRose transfer functions
//  for the HRS + 6 degree cold septum (circa 2003 - 2005).
//  R. Michaels  Sept 2008

#include "hamcTransLerColdSeptum.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTransLerColdSeptum)
#endif


hamcTransLerColdSeptum::hamcTransLerColdSeptum() 
{
}

hamcTransLerColdSeptum::~hamcTransLerColdSeptum() {
}

Int_t hamcTransLerColdSeptum::Init(hamcSpecHRS *spect) {

  which_spectrom = spect->which_spectrom;
  collim_distance = spect->GetCollimDist();

  return OK;

}

void hamcTransLerColdSeptum::Print() {
}

Int_t hamcTransLerColdSeptum::TransForm(hamcTrack *trk, Int_t where) const {

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

  Int_t origin = trk->origin;

// Functions only act from the target.

  if (origin != ITARGET) return OK;

  Int_t xbig = -9999;

  int dimen=5;
  float xtrans[5], xout[5];

  xtrans[0] = trk->tvect_orig->GetX();
  xtrans[1] = trk->tvect_orig->GetTheta();
  xtrans[2] = trk->tvect_orig->GetY();
  xtrans[3] = trk->tvect_orig->GetPhi();
  xtrans[4] = trk->tvect_orig->GetDpp();
 
  for (Int_t i=0; i<4; i++) xout[i]=xbig;
  xout[4] = xtrans[4];

// The left HRS and right HRS are slightly different

  if (which_spectrom == LEFTHRS) {

    switch(where) {

      case ISEPTIN: 

        xout[0]=x_sl_ep3_(xtrans,&dimen);
        xout[1]=t_sl_ep3_(xtrans,&dimen);
        xout[2]=y_sl_ep3_(xtrans,&dimen);
        xout[3]=p_sl_ep3_(xtrans,&dimen);
        break;

      case ISEPTOUT: 
 
        xout[0]=x_sl_ep7_(xtrans,&dimen);
        xout[1]=t_sl_ep7_(xtrans,&dimen);
        xout[2]=y_sl_ep7_(xtrans,&dimen);
        xout[3]=p_sl_ep7_(xtrans,&dimen);
        break;

      case ICOLLIM: 
      case ICOLLIM2: 
      case ICOLLIM3: 
         break;

      case IQ1EXIT: 

        xout[0]=x_sl_q1ex_(xtrans,&dimen);
        xout[1]=t_sl_q1ex_(xtrans,&dimen);
        xout[2]=y_sl_q1ex_(xtrans,&dimen);
        xout[3]=p_sl_q1ex_(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_sl_dent_(xtrans,&dimen);
        xout[1]=t_sl_dent_(xtrans,&dimen);
        xout[2]=y_sl_dent_(xtrans,&dimen);
        xout[3]=p_sl_dent_(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_sl_dext_(xtrans,&dimen);
        xout[1]=t_sl_dext_(xtrans,&dimen);
        xout[2]=y_sl_dext_(xtrans,&dimen);
        xout[3]=p_sl_dext_(xtrans,&dimen);
        break;

      case IQ3IN:

        xout[0]=x_sl_q3en_(xtrans,&dimen);
        xout[1]=t_sl_q3en_(xtrans,&dimen);
        xout[2]=y_sl_q3en_(xtrans,&dimen);
        xout[3]=p_sl_q3en_(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_sl_q3ex_(xtrans,&dimen);
        xout[1]=t_sl_q3ex_(xtrans,&dimen);
        xout[2]=y_sl_q3ex_(xtrans,&dimen);
        xout[3]=p_sl_q3ex_(xtrans,&dimen);
        break;

      case IFOCAL:
      case IPLANE1:
      case IPLANE2:

        xout[0]=x_sl_fp_(xtrans,&dimen);
        xout[1]=t_sl_fp_(xtrans,&dimen);
        xout[2]=y_sl_fp_(xtrans,&dimen);
        xout[3]=p_sl_fp_(xtrans,&dimen);
        break;

      default:

        cout << "TransLer::WARNING: undefined destination "<<where<<endl;

    }

  }

  if (which_spectrom == RIGHTHRS) {

    switch(where) {

      case ISEPTIN: 

        xout[0]=x_sr_ep3_(xtrans,&dimen);
        xout[1]=t_sr_ep3_(xtrans,&dimen);
        xout[2]=y_sr_ep3_(xtrans,&dimen);
        xout[3]=p_sr_ep3_(xtrans,&dimen);
        break;

      case ISEPTOUT: 
 
        xout[0]=x_sr_ep7_(xtrans,&dimen);
        xout[1]=t_sr_ep7_(xtrans,&dimen);
        xout[2]=y_sr_ep7_(xtrans,&dimen);
        xout[3]=p_sr_ep7_(xtrans,&dimen);
        break;

      case ICOLLIM: 
      case ICOLLIM2: 
      case ICOLLIM3: 
         break;

      case IQ1EXIT: 

        xout[0]=x_sr_q1ex_(xtrans,&dimen);
        xout[1]=t_sr_q1ex_(xtrans,&dimen);
        xout[2]=y_sr_q1ex_(xtrans,&dimen);
        xout[3]=p_sr_q1ex_(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_sr_dent_(xtrans,&dimen);
        xout[1]=t_sr_dent_(xtrans,&dimen);
        xout[2]=y_sr_dent_(xtrans,&dimen);
        xout[3]=p_sr_dent_(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_sr_dext_(xtrans,&dimen);
        xout[1]=t_sr_dext_(xtrans,&dimen);
        xout[2]=y_sr_dext_(xtrans,&dimen);
        xout[3]=p_sr_dext_(xtrans,&dimen);
        break;

      case IQ3IN:

        xout[0]=x_sr_q3en_(xtrans,&dimen);
        xout[1]=t_sr_q3en_(xtrans,&dimen);
        xout[2]=y_sr_q3en_(xtrans,&dimen);
        xout[3]=p_sr_q3en_(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_sr_q3ex_(xtrans,&dimen);
        xout[1]=t_sr_q3ex_(xtrans,&dimen);
        xout[2]=y_sr_q3ex_(xtrans,&dimen);
        xout[3]=p_sr_q3ex_(xtrans,&dimen);
        break;

      case IFOCAL:
      case IPLANE1:
      case IPLANE2:

        xout[0]=x_sr_fp_(xtrans,&dimen);
        xout[1]=t_sr_fp_(xtrans,&dimen);
        xout[2]=y_sr_fp_(xtrans,&dimen);
        xout[3]=p_sr_fp_(xtrans,&dimen);
        break;

      default:

        cout << "TransLer::WARNING: undefined destination "<<where<<endl;

    }

  }
  

  if (xout[0] != xbig) { // this means xout was loaded; if it wasn't
                         // loaded then we dont change trk->tvect.

     trk->tvect->Load(xout);

  }

  
  return OK;

}

