#ifndef ROOT_hamcTransLerWarmSeptum
#define ROOT_hamcTransLerWarmSeptum

//  hamcTransLerWarmSeptum -- LeRose's transfer functions
//  for the warm septum
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <map>

/* Here are the Fortran function prototypes */
/* This goes with prex_forward.f    */

      extern "C" float  x_sp_sen_(float *, int *);
      extern "C" float  t_sp_sen_(float *, int *);
      extern "C" float  y_sp_sen_(float *, int *);
      extern "C" float  p_sp_sen_(float *, int *);
      extern "C" float  l_sp_sen_(float *, int *);
      extern "C" float  x_sp_sm_(float *, int *);
      extern "C" float  t_sp_sm_(float *, int *);
      extern "C" float  y_sp_sm_(float *, int *);
      extern "C" float  p_sp_sm_(float *, int *);
      extern "C" float  l_sp_sm_(float *, int *);
      extern "C" float  x_sp_sex_(float *, int *);
      extern "C" float  t_sp_sex_(float *, int *);
      extern "C" float  y_sp_sex_(float *, int *);
      extern "C" float  p_sp_sex_(float *, int *);
      extern "C" float  l_sp_sex_(float *, int *);
      extern "C" float  x_sp_q1ex_(float *, int *);
      extern "C" float  t_sp_q1ex_(float *, int *);
      extern "C" float  y_sp_q1ex_(float *, int *);
      extern "C" float  p_sp_q1ex_(float *, int *);
      extern "C" float  l_sp_q1ex_(float *, int *);
      extern "C" float  x_sp_q2ex_(float *, int *);
      extern "C" float  t_sp_q2ex_(float *, int *);
      extern "C" float  y_sp_q2ex_(float *, int *);
      extern "C" float  p_sp_q2ex_(float *, int *);
      extern "C" float  l_sp_q2ex_(float *, int *);
      extern "C" float  x_sp_den_(float *, int *);
      extern "C" float  t_sp_den_(float *, int *);
      extern "C" float  y_sp_den_(float *, int *);
      extern "C" float  p_sp_den_(float *, int *);
      extern "C" float  l_sp_den_(float *, int *);
      extern "C" float  x_sp_dex_(float *, int *);
      extern "C" float  t_sp_dex_(float *, int *);
      extern "C" float  y_sp_dex_(float *, int *);
      extern "C" float  p_sp_dex_(float *, int *);
      extern "C" float  l_sp_dex_(float *, int *);
      extern "C" float  x_sp_q3en_(float *, int *);
      extern "C" float  t_sp_q3en_(float *, int *);
      extern "C" float  y_sp_q3en_(float *, int *);
      extern "C" float  p_sp_q3en_(float *, int *);
      extern "C" float  l_sp_q3en_(float *, int *);
      extern "C" float  x_sp_q3m_(float *, int *);
      extern "C" float  t_sp_q3m_(float *, int *);
      extern "C" float  y_sp_q3m_(float *, int *);
      extern "C" float  p_sp_q3m_(float *, int *);
      extern "C" float  l_sp_q3m_(float *, int *);
      extern "C" float  x_sp_q3ex_(float *, int *);
      extern "C" float  t_sp_q3ex_(float *, int *);
      extern "C" float  y_sp_q3ex_(float *, int *);
      extern "C" float  p_sp_q3ex_(float *, int *);
      extern "C" float  l_sp_q3ex_(float *, int *);
      extern "C" float  x_sp_fp_(float *, int *);
      extern "C" float  t_sp_fp_(float *, int *);
      extern "C" float  y_sp_fp_(float *, int *);
      extern "C" float  p_sp_fp_(float *, int *);
      extern "C" float  l_sp_fp_(float *, int *);
      extern "C" float  x_sp_col_(float *, int *);
      extern "C" float  t_sp_col_(float *, int *);
      extern "C" float  y_sp_col_(float *, int *);
      extern "C" float  p_sp_col_(float *, int *);
      extern "C" float  l_sp_col_(float *, int *);
      extern "C" float  x_sp_cq1x_(float *, int *);
      extern "C" float  t_sp_cq1x_(float *, int *);
      extern "C" float  y_sp_cq1x_(float *, int *);
      extern "C" float  p_sp_cq1x_(float *, int *);
      extern "C" float  l_sp_cq1x_(float *, int *);
      extern "C" float  x_sp_cden_(float *, int *);
      extern "C" float  t_sp_cden_(float *, int *);
      extern "C" float  y_sp_cden_(float *, int *);
      extern "C" float  p_sp_cden_(float *, int *);
      extern "C" float  l_sp_cden_(float *, int *);
      extern "C" float  x_sp_cdex_(float *, int *);
      extern "C" float  t_sp_cdex_(float *, int *);
      extern "C" float  y_sp_cdex_(float *, int *);
      extern "C" float  p_sp_cdex_(float *, int *);
      extern "C" float  l_sp_cdex_(float *, int *);
      extern "C" float  x_sp_cq3e_(float *, int *);
      extern "C" float  t_sp_cq3e_(float *, int *);
      extern "C" float  y_sp_cq3e_(float *, int *);
      extern "C" float  p_sp_cq3e_(float *, int *);
      extern "C" float  l_sp_cq3e_(float *, int *);
      extern "C" float  x_sp_cq3x_(float *, int *);
      extern "C" float  t_sp_cq3x_(float *, int *);
      extern "C" float  y_sp_cq3x_(float *, int *);
      extern "C" float  p_sp_cq3x_(float *, int *);
      extern "C" float  l_sp_cq3x_(float *, int *);
      extern "C" float  x_sp_cfp_(float *, int *);
      extern "C" float  t_sp_cfp_(float *, int *);
      extern "C" float  y_sp_cfp_(float *, int *);
      extern "C" float  p_sp_cfp_(float *, int *);
      extern "C" float  l_sp_cfp_(float *, int *);


class hamcSpecHRS;
 
class hamcTransLerWarmSeptum : public hamcTrans {

  public:

     hamcTransLerWarmSeptum();
     virtual ~hamcTransLerWarmSeptum(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

  protected:

  private: 

     hamcTransLerWarmSeptum(const hamcTransLerWarmSeptum& trans);
     hamcTransLerWarmSeptum& operator=(const hamcTransLerWarmSeptum& trans);

#ifndef NODICT
ClassDef (hamcTransLerWarmSeptum, 0)   // LeRose Transfer Function, Warm Septum
#endif
  };

#endif



   
