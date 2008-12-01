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

      extern "C" double  x_sp_sen__(float *, long *);
      extern "C" double  t_sp_sen__(float *, long *);
      extern "C" double  y_sp_sen__(float *, long *);
      extern "C" double  p_sp_sen__(float *, long *);
      extern "C" double  l_sp_sen__(float *, long *);
      extern "C" double  x_sp_sm__(float *, long *);
      extern "C" double  t_sp_sm__(float *, long *);
      extern "C" double  y_sp_sm__(float *, long *);
      extern "C" double  p_sp_sm__(float *, long *);
      extern "C" double  l_sp_sm__(float *, long *);
      extern "C" double  x_sp_sex__(float *, long *);
      extern "C" double  t_sp_sex__(float *, long *);
      extern "C" double  y_sp_sex__(float *, long *);
      extern "C" double  p_sp_sex__(float *, long *);
      extern "C" double  l_sp_sex__(float *, long *);
      extern "C" double  x_sp_q1ex__(float *, long *);
      extern "C" double  t_sp_q1ex__(float *, long *);
      extern "C" double  y_sp_q1ex__(float *, long *);
      extern "C" double  p_sp_q1ex__(float *, long *);
      extern "C" double  l_sp_q1ex__(float *, long *);
      extern "C" double  x_sp_q2ex__(float *, long *);
      extern "C" double  t_sp_q2ex__(float *, long *);
      extern "C" double  y_sp_q2ex__(float *, long *);
      extern "C" double  p_sp_q2ex__(float *, long *);
      extern "C" double  l_sp_q2ex__(float *, long *);
      extern "C" double  x_sp_den__(float *, long *);
      extern "C" double  t_sp_den__(float *, long *);
      extern "C" double  y_sp_den__(float *, long *);
      extern "C" double  p_sp_den__(float *, long *);
      extern "C" double  l_sp_den__(float *, long *);
      extern "C" double  x_sp_dex__(float *, long *);
      extern "C" double  t_sp_dex__(float *, long *);
      extern "C" double  y_sp_dex__(float *, long *);
      extern "C" double  p_sp_dex__(float *, long *);
      extern "C" double  l_sp_dex__(float *, long *);
      extern "C" double  x_sp_q3en__(float *, long *);
      extern "C" double  t_sp_q3en__(float *, long *);
      extern "C" double  y_sp_q3en__(float *, long *);
      extern "C" double  p_sp_q3en__(float *, long *);
      extern "C" double  l_sp_q3en__(float *, long *);
      extern "C" double  x_sp_q3m__(float *, long *);
      extern "C" double  t_sp_q3m__(float *, long *);
      extern "C" double  y_sp_q3m__(float *, long *);
      extern "C" double  p_sp_q3m__(float *, long *);
      extern "C" double  l_sp_q3m__(float *, long *);
      extern "C" double  x_sp_q3ex__(float *, long *);
      extern "C" double  t_sp_q3ex__(float *, long *);
      extern "C" double  y_sp_q3ex__(float *, long *);
      extern "C" double  p_sp_q3ex__(float *, long *);
      extern "C" double  l_sp_q3ex__(float *, long *);
      extern "C" double  x_sp_fp__(float *, long *);
      extern "C" double  t_sp_fp__(float *, long *);
      extern "C" double  y_sp_fp__(float *, long *);
      extern "C" double  p_sp_fp__(float *, long *);
      extern "C" double  l_sp_fp__(float *, long *);
      extern "C" double  x_sp_col__(float *, long *);
      extern "C" double  t_sp_col__(float *, long *);
      extern "C" double  y_sp_col__(float *, long *);
      extern "C" double  p_sp_col__(float *, long *);
      extern "C" double  l_sp_col__(float *, long *);
      extern "C" double  x_sp_cq1x__(float *, long *);
      extern "C" double  t_sp_cq1x__(float *, long *);
      extern "C" double  y_sp_cq1x__(float *, long *);
      extern "C" double  p_sp_cq1x__(float *, long *);
      extern "C" double  l_sp_cq1x__(float *, long *);
      extern "C" double  x_sp_cden__(float *, long *);
      extern "C" double  t_sp_cden__(float *, long *);
      extern "C" double  y_sp_cden__(float *, long *);
      extern "C" double  p_sp_cden__(float *, long *);
      extern "C" double  l_sp_cden__(float *, long *);
      extern "C" double  x_sp_cdex__(float *, long *);
      extern "C" double  t_sp_cdex__(float *, long *);
      extern "C" double  y_sp_cdex__(float *, long *);
      extern "C" double  p_sp_cdex__(float *, long *);
      extern "C" double  l_sp_cdex__(float *, long *);
      extern "C" double  x_sp_cq3e__(float *, long *);
      extern "C" double  t_sp_cq3e__(float *, long *);
      extern "C" double  y_sp_cq3e__(float *, long *);
      extern "C" double  p_sp_cq3e__(float *, long *);
      extern "C" double  l_sp_cq3e__(float *, long *);
      extern "C" double  x_sp_cq3x__(float *, long *);
      extern "C" double  t_sp_cq3x__(float *, long *);
      extern "C" double  y_sp_cq3x__(float *, long *);
      extern "C" double  p_sp_cq3x__(float *, long *);
      extern "C" double  l_sp_cq3x__(float *, long *);
      extern "C" double  x_sp_cfp__(float *, long *);
      extern "C" double  t_sp_cfp__(float *, long *);
      extern "C" double  y_sp_cfp__(float *, long *);
      extern "C" double  p_sp_cfp__(float *, long *);
      extern "C" double  l_sp_cfp__(float *, long *);


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



   
