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

      extern "C" double  x_s4_sext__(float *, long *);
      extern "C" double  t_s4_sext__(float *, long *);
      extern "C" double  y_s4_sext__(float *, long *);
      extern "C" double  p_s4_sext__(float *, long *);
      extern "C" double  l_s4_sext__(float *, long *);
      extern "C" double  x_s4_q1en__(float *, long *);
      extern "C" double  t_s4_q1en__(float *, long *);
      extern "C" double  y_s4_q1en__(float *, long *);
      extern "C" double  p_s4_q1en__(float *, long *);
      extern "C" double  l_s4_q1en__(float *, long *);
      extern "C" double  x_s4_q1ex__(float *, long *);
      extern "C" double  t_s4_q1ex__(float *, long *);
      extern "C" double  y_s4_q1ex__(float *, long *);
      extern "C" double  p_s4_q1ex__(float *, long *);
      extern "C" double  l_s4_q1ex__(float *, long *);
      extern "C" double  x_s4_den__(float *, long *);
      extern "C" double  t_s4_den__(float *, long *);
      extern "C" double  y_s4_den__(float *, long *);
      extern "C" double  p_s4_den__(float *, long *);
      extern "C" double  l_s4_den__(float *, long *);
      extern "C" double  x_s4_dex__(float *, long *);
      extern "C" double  t_s4_dex__(float *, long *);
      extern "C" double  y_s4_dex__(float *, long *);
      extern "C" double  p_s4_dex__(float *, long *);
      extern "C" double  l_s4_dex__(float *, long *);
      extern "C" double  x_s4_q3en__(float *, long *);
      extern "C" double  t_s4_q3en__(float *, long *);
      extern "C" double  y_s4_q3en__(float *, long *);
      extern "C" double  p_s4_q3en__(float *, long *);
      extern "C" double  l_s4_q3en__(float *, long *);
      extern "C" double  x_s4_q3ex__(float *, long *);
      extern "C" double  t_s4_q3ex__(float *, long *);
      extern "C" double  y_s4_q3ex__(float *, long *);
      extern "C" double  p_s4_q3ex__(float *, long *);
      extern "C" double  l_s4_q3ex__(float *, long *);
      extern "C" double  x_s4_fp__(float *, long *);
      extern "C" double  t_s4_fp__(float *, long *);
      extern "C" double  y_s4_fp__(float *, long *);
      extern "C" double  p_s4_fp__(float *, long *);
      extern "C" double  l_s4_fp__(float *, long *);

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



   
