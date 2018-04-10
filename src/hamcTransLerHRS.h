#ifndef ROOT_hamcTransLerHRS
#define ROOT_hamcTransLerHRS

//  hamcTransLerHRS -- LeRose's transfer functions
//  for the High Resolution Spectrometers (standard setup) 
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <map>

// Functions for standard HRS.  from monte_trans_hrs.f
// "e" means left HRS, "h" means right HRS.

      extern "C" float  x_e_q1ex_(float *, int *);
      extern "C" float  t_e_q1ex_(float *, int *);
      extern "C" float  y_e_q1ex_(float *, int *);
      extern "C" float  p_e_q1ex_(float *, int *);
      extern "C" float  l_e_q1ex_(float *, int *);
      extern "C" float  x_e_dent_(float *, int *);
      extern "C" float  t_e_dent_(float *, int *);
      extern "C" float  y_e_dent_(float *, int *);
      extern "C" float  p_e_dent_(float *, int *);
      extern "C" float  l_e_dent_(float *, int *);
      extern "C" float  x_e_dext_(float *, int *);
      extern "C" float  t_e_dext_(float *, int *);
      extern "C" float  y_e_dext_(float *, int *);
      extern "C" float  p_e_dext_(float *, int *);
      extern "C" float  l_e_dext_(float *, int *);
      extern "C" float  x_e_q3en_(float *, int *);
      extern "C" float  t_e_q3en_(float *, int *);
      extern "C" float  y_e_q3en_(float *, int *);
      extern "C" float  p_e_q3en_(float *, int *);
      extern "C" float  l_e_q3en_(float *, int *);
      extern "C" float  x_e_q3ex_(float *, int *);
      extern "C" float  t_e_q3ex_(float *, int *);
      extern "C" float  y_e_q3ex_(float *, int *);
      extern "C" float  p_e_q3ex_(float *, int *);
      extern "C" float  l_e_q3ex_(float *, int *);
      extern "C" float  x_e_fp_(float *, int *);
      extern "C" float  t_e_fp_(float *, int *);
      extern "C" float  y_e_fp_(float *, int *);
      extern "C" float  p_e_fp_(float *, int *);
      extern "C" float  l_e_fp_(float *, int *);
      extern "C" float  x_h_q1ex_(float *, int *);
      extern "C" float  t_h_q1ex_(float *, int *);
      extern "C" float  y_h_q1ex_(float *, int *);
      extern "C" float  p_h_q1ex_(float *, int *);
      extern "C" float  l_h_q1ex_(float *, int *);
      extern "C" float  x_h_dent_(float *, int *);
      extern "C" float  t_h_dent_(float *, int *);
      extern "C" float  y_h_dent_(float *, int *);
      extern "C" float  p_h_dent_(float *, int *);
      extern "C" float  l_h_dent_(float *, int *);
      extern "C" float  x_h_dext_(float *, int *);
      extern "C" float  t_h_dext_(float *, int *);
      extern "C" float  y_h_dext_(float *, int *);
      extern "C" float  p_h_dext_(float *, int *);
      extern "C" float  l_h_dext_(float *, int *);
      extern "C" float  x_h_q3en_(float *, int *);
      extern "C" float  t_h_q3en_(float *, int *);
      extern "C" float  y_h_q3en_(float *, int *);
      extern "C" float  p_h_q3en_(float *, int *);
      extern "C" float  l_h_q3en_(float *, int *);
      extern "C" float  x_h_q3ex_(float *, int *);
      extern "C" float  t_h_q3ex_(float *, int *);
      extern "C" float  y_h_q3ex_(float *, int *);
      extern "C" float  p_h_q3ex_(float *, int *);
      extern "C" float  l_h_q3ex_(float *, int *);
      extern "C" float  x_h_fp_(float *, int *);
      extern "C" float  t_h_fp_(float *, int *);
      extern "C" float  y_h_fp_(float *, int *);
      extern "C" float  p_h_fp_(float *, int *);
      extern "C" float  l_h_fp_(float *, int *);


class hamcSpecHRS;
 
class hamcTransLerHRS : public hamcTrans {

  public:

     hamcTransLerHRS();
     virtual ~hamcTransLerHRS(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

  protected:

  private: 

     hamcTransLerHRS(const hamcTransLerHRS& trans);
     hamcTransLerHRS& operator=(const hamcTransLerHRS& trans);

#ifndef NODICT
ClassDef (hamcTransLerHRS, 0)   // LeRose Transfer Function, std. HRS
#endif
  };

#endif



   
