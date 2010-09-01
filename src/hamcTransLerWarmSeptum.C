//  hamcTransLerWarmSeptum   -- Transport model using LeRose transfer functions
//  for the new warm 5 degree septum for PREX 
//  R. Michaels  Sept 2008

#include "hamcTransLerWarmSeptum.h"
#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTransLerWarmSeptum)
#endif


hamcTransLerWarmSeptum::hamcTransLerWarmSeptum() 
{
}

hamcTransLerWarmSeptum::~hamcTransLerWarmSeptum() {
}

Int_t hamcTransLerWarmSeptum::Init(hamcSpecHRS *spect) {

  which_spectrom = spect->which_spectrom;
  collim_distance = spect->GetCollimDist();

  return OK;

}

void hamcTransLerWarmSeptum::Print() {
}

Int_t hamcTransLerWarmSeptum::TransForm(hamcTrack *trk, Int_t where) const {

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
  Int_t debug_trans = 0; // to get do some debugging here.

  long dimen=5;
  float xtrans[5], xout[5];

  xtrans[0] = trk->tvect_orig->GetX();
  xtrans[1] = trk->tvect_orig->GetTheta();
  xtrans[2] = trk->tvect_orig->GetY();
  xtrans[3] = trk->tvect_orig->GetPhi();
  xtrans[4] = trk->tvect_orig->GetDpp();

  for (Int_t i=0; i<4; i++) xout[i]=xbig;
  xout[4] = xtrans[4];

  if (debug_trans &&  (origin == ITARGET && xtrans[1] > -0.04 && xtrans[1] < 0.04 && xtrans[3] > -0.04 && xtrans[3] < 0.04) && xtrans[4] > -0.02 && xtrans[4] < 0.02) {
// Test: compare (tgt -> fp) to (tgt -> col -> fp) 
// idiot test:     xtrans[0] = 0.;
//     xtrans[1] = 0.02;
//     xtrans[2] = 0.;
//     xtrans[3] = 0.02;
//     xtrans[4] = 0.001;
     Float_t xtof=x_sp_fp__(xtrans,&dimen);    // tgt -> fp
     Float_t ttof=t_sp_fp__(xtrans,&dimen);
     Float_t ytof=y_sp_fp__(xtrans,&dimen);
     Float_t ptof=p_sp_fp__(xtrans,&dimen);     
     xout[0]=x_sp_col__(xtrans,&dimen);// tgt -> col
     xout[1]=t_sp_col__(xtrans,&dimen);
     xout[2]=y_sp_col__(xtrans,&dimen);
     xout[3]=p_sp_col__(xtrans,&dimen);
     xout[4]=xtrans[4];
     Float_t xcolf=x_sp_cfp__(xout,&dimen);    // col -> fp
     Float_t tcolf=t_sp_cfp__(xout,&dimen);
     Float_t ycolf=y_sp_cfp__(xout,&dimen);
     Float_t pcolf=p_sp_cfp__(xout,&dimen);
     cout << endl << endl;
     cout <<"orig flag "<<origin<<endl;
     cout <<"orig "<<xtrans[0]<<" "<<xtrans[1]<<" "<<xtrans[2]<<" "<<xtrans[3]<<" "<<xtrans[4]<<endl;
     cout <<"T->F "<<xtof<<"  "<<ttof<<"  "<<ytof<<"  "<<ptof<<endl;
     cout <<"T->C "<<xcolf<<"  "<<tcolf<<"  "<<ycolf<<"  "<<pcolf<<endl;
  }

  if (origin == ITARGET) {

    switch(where) {

      case ISEPTIN: 

        xout[0]=x_sp_sen__(xtrans,&dimen);
        xout[1]=t_sp_sen__(xtrans,&dimen);
        xout[2]=y_sp_sen__(xtrans,&dimen);
        xout[3]=p_sp_sen__(xtrans,&dimen);
        break;

      case ISEPTOUT: 
 
        xout[0]=x_sp_sex__(xtrans,&dimen);
        xout[1]=t_sp_sex__(xtrans,&dimen);
        xout[2]=y_sp_sex__(xtrans,&dimen);
        xout[3]=p_sp_sex__(xtrans,&dimen);
        break;

      case ICOLLIM3: 

	trk->tvect->Load(xtrans);
        Drift(0.797, trk);  // 79.7 cm to sieve.
        return OK;

      case ICOLLIM: 
      case ICOLLIM2: 

        xout[0]=x_sp_col__(xtrans,&dimen);
        xout[1]=t_sp_col__(xtrans,&dimen);
        xout[2]=y_sp_col__(xtrans,&dimen);
        xout[3]=p_sp_col__(xtrans,&dimen);
        break;

      case IQ1EXIT: 

        xout[0]=x_sp_q1ex__(xtrans,&dimen);
        xout[1]=t_sp_q1ex__(xtrans,&dimen);
        xout[2]=y_sp_q1ex__(xtrans,&dimen);
        xout[3]=p_sp_q1ex__(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_sp_den__(xtrans,&dimen);
        xout[1]=t_sp_den__(xtrans,&dimen);
        xout[2]=y_sp_den__(xtrans,&dimen);
        xout[3]=p_sp_den__(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_sp_dex__(xtrans,&dimen);
        xout[1]=t_sp_dex__(xtrans,&dimen);
        xout[2]=y_sp_dex__(xtrans,&dimen);
        xout[3]=p_sp_dex__(xtrans,&dimen);

        break;

      case IQ3IN:

        xout[0]=x_sp_q3en__(xtrans,&dimen);
        xout[1]=t_sp_q3en__(xtrans,&dimen);
        xout[2]=y_sp_q3en__(xtrans,&dimen);
        xout[3]=p_sp_q3en__(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_sp_q3ex__(xtrans,&dimen);
        xout[1]=t_sp_q3ex__(xtrans,&dimen);
        xout[2]=y_sp_q3ex__(xtrans,&dimen);
        xout[3]=p_sp_q3ex__(xtrans,&dimen);
        break;

      case IFOCAL:

        xout[0]=x_sp_fp__(xtrans,&dimen);
        xout[1]=t_sp_fp__(xtrans,&dimen);
        xout[2]=y_sp_fp__(xtrans,&dimen);
        xout[3]=p_sp_fp__(xtrans,&dimen);
        break;

      case IPLANE1:

	xout[0]=x_sp_fp__(xtrans,&dimen);
        xout[1]=t_sp_fp__(xtrans,&dimen);
        xout[2]=y_sp_fp__(xtrans,&dimen);
        xout[3]=p_sp_fp__(xtrans,&dimen);

        trk->tvect->Load(xout);
        Drift(1.0, trk);  // 1 meter
        return OK;

      case IPLANE2:

	xout[0]=x_sp_fp__(xtrans,&dimen);
        xout[1]=t_sp_fp__(xtrans,&dimen);
        xout[2]=y_sp_fp__(xtrans,&dimen);
        xout[3]=p_sp_fp__(xtrans,&dimen);

        trk->tvect->Load(xout);
        Drift(1.48, trk);  // 1.48 m
        return OK;

      default:

        cout << "TransLer::WARNING: undefined destination "<<where<<endl;

    }

  }
  
  if (origin == ICOLLIM || origin == ICOLLIM2 || origin == ICOLLIM3) {

    switch(where) {

      case IQ1EXIT: 

        xout[0]=x_sp_cq1x__(xtrans,&dimen);
        xout[1]=t_sp_cq1x__(xtrans,&dimen);
        xout[2]=y_sp_cq1x__(xtrans,&dimen);
        xout[3]=p_sp_cq1x__(xtrans,&dimen);
        break;

      case IDIPIN:

        xout[0]=x_sp_cden__(xtrans,&dimen);
        xout[1]=t_sp_cden__(xtrans,&dimen);
        xout[2]=y_sp_cden__(xtrans,&dimen);
        xout[3]=p_sp_cden__(xtrans,&dimen);
        break;

      case IDIPEXIT:

        xout[0]=x_sp_cdex__(xtrans,&dimen);
        xout[1]=t_sp_cdex__(xtrans,&dimen);
        xout[2]=y_sp_cdex__(xtrans,&dimen);
        xout[3]=p_sp_cdex__(xtrans,&dimen);
        break;

      case IQ3IN:

        xout[0]=x_sp_cq3e__(xtrans,&dimen);
        xout[1]=t_sp_cq3e__(xtrans,&dimen);
        xout[2]=y_sp_cq3e__(xtrans,&dimen);
        xout[3]=p_sp_cq3e__(xtrans,&dimen);
        break;

      case IQ3EXIT:

        xout[0]=x_sp_cq3x__(xtrans,&dimen);
        xout[1]=t_sp_cq3x__(xtrans,&dimen);
        xout[2]=y_sp_cq3x__(xtrans,&dimen);
        xout[3]=p_sp_cq3x__(xtrans,&dimen);
        break;

      case IFOCAL:

        xout[0]=x_sp_cfp__(xtrans,&dimen);
        xout[1]=t_sp_cfp__(xtrans,&dimen);
        xout[2]=y_sp_cfp__(xtrans,&dimen);
        xout[3]=p_sp_cfp__(xtrans,&dimen);
        break;

      case IPLANE1:

	xout[0]=x_sp_cfp__(xtrans,&dimen);
        xout[1]=t_sp_cfp__(xtrans,&dimen);
        xout[2]=y_sp_cfp__(xtrans,&dimen);
        xout[3]=p_sp_cfp__(xtrans,&dimen);

        trk->tvect->Load(xout);
        Drift(1.0, trk);  // 1 meter
        return OK;

      case IPLANE2:

	xout[0]=x_sp_cfp__(xtrans,&dimen);
        xout[1]=t_sp_cfp__(xtrans,&dimen);
        xout[2]=y_sp_cfp__(xtrans,&dimen);
        xout[3]=p_sp_cfp__(xtrans,&dimen);

        trk->tvect->Load(xout);
        Drift(1.48, trk);  // 1.48 m
        return OK;

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

