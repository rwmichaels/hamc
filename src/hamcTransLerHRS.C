//  hamcTransLerHRS   -- Transport model using LeRose transfer functions
//  for the standard HRS setup.
//  R. Michaels  Sept 2008

#include "hamcTransLerHRS.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTransLerHRS)
#endif


hamcTransLerHRS::hamcTransLerHRS() 
{
  collim_distance=0;
}

hamcTransLerHRS::~hamcTransLerHRS() {
}

Int_t hamcTransLerHRS::Init(hamcSpecHRS *spect) {

  which_spectrom = spect->which_spectrom;
  collim_distance = spect->GetCollimDist();
  return OK;

}

void hamcTransLerHRS::Print() {
}

Int_t hamcTransLerHRS::TransForm(hamcTrack *trk, Int_t where) const {

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
      case ISEPTOUT: 
	
        return OK;    // There is no septum.

      case ICOLLIM: 
      case ICOLLIM2: 
      case ICOLLIM3: 
 
        Drift(collim_distance, trk);
        return OK;

      case IQ1EXIT: 

        xout[0]=x_e_q1ex_(xtrans,&dimen);
        xout[1]=t_e_q1ex_(xtrans,&dimen);
        xout[2]=y_e_q1ex_(xtrans,&dimen);
        xout[3]=p_e_q1ex_(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_e_dent_(xtrans,&dimen);
        xout[1]=t_e_dent_(xtrans,&dimen);
        xout[2]=y_e_dent_(xtrans,&dimen);
        xout[3]=p_e_dent_(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_e_dext_(xtrans,&dimen);
        xout[1]=t_e_dext_(xtrans,&dimen);
        xout[2]=y_e_dext_(xtrans,&dimen);
        xout[3]=p_e_dext_(xtrans,&dimen);
        break;

      case IQ3IN:

        xout[0]=x_e_q3en_(xtrans,&dimen);
        xout[1]=t_e_q3en_(xtrans,&dimen);
        xout[2]=y_e_q3en_(xtrans,&dimen);
        xout[3]=p_e_q3en_(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_e_q3ex_(xtrans,&dimen);
        xout[1]=t_e_q3ex_(xtrans,&dimen);
        xout[2]=y_e_q3ex_(xtrans,&dimen);
        xout[3]=p_e_q3ex_(xtrans,&dimen);
        break;

      case IFOCAL:
      case IPLANE1:
      case IPLANE2:

        xout[0]=x_e_fp_(xtrans,&dimen);
        xout[1]=t_e_fp_(xtrans,&dimen);
        xout[2]=y_e_fp_(xtrans,&dimen);
        xout[3]=p_e_fp_(xtrans,&dimen);
        break;

      default:

        cout << "TransLer::WARNING: undefined destination"<<endl;

    }

  }

  if (which_spectrom == RIGHTHRS) {

    switch(where) {

      case ISEPTIN: 
      case ISEPTOUT: 
	
        return OK;    // There is no septum.

      case ICOLLIM: 
      case ICOLLIM2: 
      case ICOLLIM3: 
 
        Drift(collim_distance, trk);
        return OK;

      case IQ1EXIT: 

        xout[0]=x_h_q1ex_(xtrans,&dimen);
        xout[1]=t_h_q1ex_(xtrans,&dimen);
        xout[2]=y_h_q1ex_(xtrans,&dimen);
        xout[3]=p_h_q1ex_(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_h_dent_(xtrans,&dimen);
        xout[1]=t_h_dent_(xtrans,&dimen);
        xout[2]=y_h_dent_(xtrans,&dimen);
        xout[3]=p_h_dent_(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_h_dext_(xtrans,&dimen);
        xout[1]=t_h_dext_(xtrans,&dimen);
        xout[2]=y_h_dext_(xtrans,&dimen);
        xout[3]=p_h_dext_(xtrans,&dimen);
        break;

      case IQ3IN:

        xout[0]=x_h_q3en_(xtrans,&dimen);
        xout[1]=t_h_q3en_(xtrans,&dimen);
        xout[2]=y_h_q3en_(xtrans,&dimen);
        xout[3]=p_h_q3en_(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_h_q3ex_(xtrans,&dimen);
        xout[1]=t_h_q3ex_(xtrans,&dimen);
        xout[2]=y_h_q3ex_(xtrans,&dimen);
        xout[3]=p_h_q3ex_(xtrans,&dimen);
        break;

      case IFOCAL:
      case IPLANE1:
      case IPLANE2:

        xout[0]=x_h_fp_(xtrans,&dimen);
        xout[1]=t_h_fp_(xtrans,&dimen);
        xout[2]=y_h_fp_(xtrans,&dimen);
        xout[3]=p_h_fp_(xtrans,&dimen);
        break;

      default:

        cout << "TransLer::WARNING: undefined destination"<<endl;

    }

  }

  //  cout << "\n---- hamcTransLerHRS "<<where<<endl;
  //  cout << "    in "<<xtrans[0]<<"  "<<xtrans[1]<<"  "<<xtrans[2]<<"  "<<xtrans[3]<<endl;
  //  cout << "    out "<<xout[0]<<"  "<<xout[1]<<"  "<<xout[2]<<"  "<<xout[3]<<endl;

  if (xout[0] != xbig) { // this means xout was loaded; if it wasn't
                         // loaded then we dont change trk->tvect.

     trk->tvect->Load(xout);

  }

  return OK;

}

