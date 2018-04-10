#ifndef ROOT_hamcTransLerColdSeptum
#define ROOT_hamcTransLerColdSeptum

//  hamcTransLerColdSeptum -- LeRose's transfer functions
//  for the cold septum (circa 2003-2005).
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <map>

// Functions for HRS + 6 degree septum.  
// sl means left HRS,  sr means Right HRS  (s means septum)

// First, the Left HRS from ls_6d_forward.f

      extern "C" float  x_sl_ep3_(float *, int *);
      extern "C" float  t_sl_ep3_(float *, int *);
      extern "C" float  y_sl_ep3_(float *, int *);
      extern "C" float  p_sl_ep3_(float *, int *);
      extern "C" float  l_sl_ep3_(float *, int *);
      extern "C" float  x_sl_ep4_(float *, int *);
      extern "C" float  t_sl_ep4_(float *, int *);
      extern "C" float  y_sl_ep4_(float *, int *);
      extern "C" float  p_sl_ep4_(float *, int *);
      extern "C" float  l_sl_ep4_(float *, int *);
      extern "C" float  x_sl_ep5_(float *, int *);
      extern "C" float  t_sl_ep5_(float *, int *);
      extern "C" float  y_sl_ep5_(float *, int *);
      extern "C" float  p_sl_ep5_(float *, int *);
      extern "C" float  l_sl_ep5_(float *, int *);
      extern "C" float  x_sl_ep6_(float *, int *);
      extern "C" float  t_sl_ep6_(float *, int *);
      extern "C" float  y_sl_ep6_(float *, int *);
      extern "C" float  p_sl_ep6_(float *, int *);
      extern "C" float  l_sl_ep6_(float *, int *);
      extern "C" float  x_sl_ep7_(float *, int *);
      extern "C" float  t_sl_ep7_(float *, int *);
      extern "C" float  y_sl_ep7_(float *, int *);
      extern "C" float  p_sl_ep7_(float *, int *);
      extern "C" float  l_sl_ep7_(float *, int *);
      extern "C" float  x_sl_q1ex_(float *, int *);  
      extern "C" float  t_sl_q1ex_(float *, int *);
      extern "C" float  y_sl_q1ex_(float *, int *);
      extern "C" float  p_sl_q1ex_(float *, int *);
      extern "C" float  l_sl_q1ex_(float *, int *);
      extern "C" float  x_sl_dent_(float *, int *);
      extern "C" float  t_sl_dent_(float *, int *);
      extern "C" float  y_sl_dent_(float *, int *);
      extern "C" float  p_sl_dent_(float *, int *);
      extern "C" float  l_sl_dent_(float *, int *);
      extern "C" float  x_sl_dext_(float *, int *);
      extern "C" float  t_sl_dext_(float *, int *);
      extern "C" float  y_sl_dext_(float *, int *);
      extern "C" float  p_sl_dext_(float *, int *);
      extern "C" float  l_sl_dext_(float *, int *);
      extern "C" float  x_sl_q3en_(float *, int *);
      extern "C" float  t_sl_q3en_(float *, int *);
      extern "C" float  y_sl_q3en_(float *, int *);
      extern "C" float  p_sl_q3en_(float *, int *);
      extern "C" float  l_sl_q3en_(float *, int *);
      extern "C" float  x_sl_q3ex_(float *, int *);
      extern "C" float  t_sl_q3ex_(float *, int *);
      extern "C" float  y_sl_q3ex_(float *, int *);
      extern "C" float  p_sl_q3ex_(float *, int *);
      extern "C" float  l_sl_q3ex_(float *, int *);
      extern "C" float  x_sl_fp_(float *, int *);
      extern "C" float  t_sl_fp_(float *, int *);
      extern "C" float  y_sl_fp_(float *, int *);
      extern "C" float  p_sl_fp_(float *, int *);
      extern "C" float  l_sl_fp_(float *, int *);
// 
// begin right-HRS + 6 degree septum
// from R6_forward.f
//
      extern "C" float  x_sr_ep3_(float *, int *);
      extern "C" float  t_sr_ep3_(float *, int *);
      extern "C" float  y_sr_ep3_(float *, int *);
      extern "C" float  p_sr_ep3_(float *, int *);
      extern "C" float  l_sr_ep3_(float *, int *);
      extern "C" float  x_sr_ep4_(float *, int *);
      extern "C" float  t_sr_ep4_(float *, int *);
      extern "C" float  y_sr_ep4_(float *, int *);
      extern "C" float  p_sr_ep4_(float *, int *);
      extern "C" float  l_sr_ep4_(float *, int *);
      extern "C" float  x_sr_ep5_(float *, int *);
      extern "C" float  t_sr_ep5_(float *, int *);
      extern "C" float  y_sr_ep5_(float *, int *);
      extern "C" float  p_sr_ep5_(float *, int *);
      extern "C" float  l_sr_ep5_(float *, int *);
      extern "C" float  x_sr_ep6_(float *, int *);
      extern "C" float  t_sr_ep6_(float *, int *);
      extern "C" float  y_sr_ep6_(float *, int *);
      extern "C" float  p_sr_ep6_(float *, int *);
      extern "C" float  l_sr_ep6_(float *, int *);
      extern "C" float  x_sr_ep7_(float *, int *);
      extern "C" float  t_sr_ep7_(float *, int *);
      extern "C" float  y_sr_ep7_(float *, int *);
      extern "C" float  p_sr_ep7_(float *, int *);
      extern "C" float  l_sr_ep7_(float *, int *);
      extern "C" float  x_sr_q1ex_(float *, int *);
      extern "C" float  t_sr_q1ex_(float *, int *);
      extern "C" float  y_sr_q1ex_(float *, int *);
      extern "C" float  p_sr_q1ex_(float *, int *);
      extern "C" float  l_sr_q1ex_(float *, int *);
      extern "C" float  x_sr_dent_(float *, int *);
      extern "C" float  t_sr_dent_(float *, int *);
      extern "C" float  y_sr_dent_(float *, int *);
      extern "C" float  p_sr_dent_(float *, int *);
      extern "C" float  l_sr_dent_(float *, int *);
      extern "C" float  x_sr_dext_(float *, int *);
      extern "C" float  t_sr_dext_(float *, int *);
      extern "C" float  y_sr_dext_(float *, int *);
      extern "C" float  p_sr_dext_(float *, int *);
      extern "C" float  l_sr_dext_(float *, int *);
      extern "C" float  x_sr_q3en_(float *, int *);
      extern "C" float  t_sr_q3en_(float *, int *);
      extern "C" float  y_sr_q3en_(float *, int *);
      extern "C" float  p_sr_q3en_(float *, int *);
      extern "C" float  l_sr_q3en_(float *, int *);
      extern "C" float  x_sr_q3ex_(float *, int *);
      extern "C" float  t_sr_q3ex_(float *, int *);
      extern "C" float  y_sr_q3ex_(float *, int *);
      extern "C" float  p_sr_q3ex_(float *, int *);
      extern "C" float  l_sr_q3ex_(float *, int *);
      extern "C" float  x_sr_fp_(float *, int *);
      extern "C" float  t_sr_fp_(float *, int *);
      extern "C" float  y_sr_fp_(float *, int *);
      extern "C" float  p_sr_fp_(float *, int *);
      extern "C" float  l_sr_fp_(float *, int *);


class hamcSpecHRS;
 
class hamcTransLerColdSeptum : public hamcTrans {

  public:

     hamcTransLerColdSeptum();
     virtual ~hamcTransLerColdSeptum(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

  protected:

  private: 

     hamcTransLerColdSeptum(const hamcTransLerColdSeptum& trans);
     hamcTransLerColdSeptum& operator=(const hamcTransLerColdSeptum& trans);

#ifndef NODICT
ClassDef (hamcTransLerColdSeptum, 0)   // LeRose Transfer Function, Cold Septum
#endif
  };

#endif



   
