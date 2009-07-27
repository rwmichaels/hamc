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

   vector<Float_t> din;
   Float_t dout[6];

   din = trk->tvect->Get();

   dout[0] = din[0] + dist*din[1];
   dout[1] = din[1];
   dout[2] = din[2] + dist*din[3];
   dout[4] = din[4];
   dout[5] = din[5];

   trk->tvect->Load(dout);

   return OK;
}
