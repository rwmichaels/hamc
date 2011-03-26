//  hamcTrans   -- Transport (abstract interface)
//  R. Michaels  June 2008

#include "hamcTrans.h"
#include "hamcTrack.h"
#include "hamcSpecHRS.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcTrans)
#endif


hamcTrans::hamcTrans() : did_init(kFALSE),which_spectrom(0)
{ }

hamcTrans::~hamcTrans() {
}


Int_t hamcTrans::Drift(Float_t dist, hamcTrack *trk) const {
// transform through a drift space 'dist'

   Float_t dout[6];

   dout[0] = trk->tvect->GetX() + dist*trk->tvect->GetTheta();
   dout[1] = trk->tvect->GetTheta();
   dout[2] = trk->tvect->GetY() + dist*trk->tvect->GetPhi();
   dout[3] = trk->tvect->GetPhi();
   dout[4] = trk->tvect->GetDpp();

   trk->tvect->Load(dout);

   return OK;
}
