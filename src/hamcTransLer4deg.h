#ifndef ROOT_hamcTransLer4deg
#define ROOT_hamcTransLer4deg

//  hamcTransLer4deg -- LeRose's transfer functions
//  for the 4 degree septum
//  R. Michaels  April 2013

#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <map>

/* Here are the Fortran function prototypes */
/* This goes with crex_4degr.f    */

      extern "C" float  x_s4_sext_(float *, int *);
      extern "C" float  t_s4_sext_(float *, int *);
      extern "C" float  y_s4_sext_(float *, int *);
      extern "C" float  p_s4_sext_(float *, int *);
      extern "C" float  l_s4_sext_(float *, int *);
      extern "C" float  x_s4_q1en_(float *, int *);
      extern "C" float  t_s4_q1en_(float *, int *);
      extern "C" float  y_s4_q1en_(float *, int *);
      extern "C" float  p_s4_q1en_(float *, int *);
      extern "C" float  l_s4_q1en_(float *, int *);
      extern "C" float  x_s4_q1ex_(float *, int *);
      extern "C" float  t_s4_q1ex_(float *, int *);
      extern "C" float  y_s4_q1ex_(float *, int *);
      extern "C" float  p_s4_q1ex_(float *, int *);
      extern "C" float  l_s4_q1ex_(float *, int *);
      extern "C" float  x_s4_den_(float *, int *);
      extern "C" float  t_s4_den_(float *, int *);
      extern "C" float  y_s4_den_(float *, int *);
      extern "C" float  p_s4_den_(float *, int *);
      extern "C" float  l_s4_den_(float *, int *);
      extern "C" float  x_s4_dex_(float *, int *);
      extern "C" float  t_s4_dex_(float *, int *);
      extern "C" float  y_s4_dex_(float *, int *);
      extern "C" float  p_s4_dex_(float *, int *);
      extern "C" float  l_s4_dex_(float *, int *);
      extern "C" float  x_s4_q3en_(float *, int *);
      extern "C" float  t_s4_q3en_(float *, int *);
      extern "C" float  y_s4_q3en_(float *, int *);
      extern "C" float  p_s4_q3en_(float *, int *);
      extern "C" float  l_s4_q3en_(float *, int *);
      extern "C" float  x_s4_q3ex_(float *, int *);
      extern "C" float  t_s4_q3ex_(float *, int *);
      extern "C" float  y_s4_q3ex_(float *, int *);
      extern "C" float  p_s4_q3ex_(float *, int *);
      extern "C" float  l_s4_q3ex_(float *, int *);
      extern "C" float  x_s4_fp_(float *, int *);
      extern "C" float  t_s4_fp_(float *, int *);
      extern "C" float  y_s4_fp_(float *, int *);
      extern "C" float  p_s4_fp_(float *, int *);
      extern "C" float  l_s4_fp_(float *, int *);

class hamcSpecHRS;
 
class hamcTransLer4deg : public hamcTrans {

  public:

     hamcTransLer4deg();
     virtual ~hamcTransLer4deg(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

  protected:

  private: 

     hamcTransLer4deg(const hamcTransLer4deg& trans);
     hamcTransLer4deg& operator=(const hamcTransLer4deg& trans);

#ifndef NODICT
ClassDef (hamcTransLer4deg, 0)   // LeRose Transfer Function, 4 degr. Septum
#endif
  };

#endif



   
