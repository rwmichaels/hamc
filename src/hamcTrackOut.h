#ifndef ROOT_hamcTrackOut
#define ROOT_hamcTrackOut

//  hamcTrackOut   -- class for a track going out from target
//  R. Michaels  Apr 2008

#include "Rtypes.h"
#include "hamcTrack.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include <vector>
#include <map>

class hamcBeam;

class hamcTrackOut : public hamcTrack {

  public:

     hamcTrackOut();
     hamcTrackOut(std::string pid, Float_t energy, Float_t x, Float_t theta, Float_t y, Float_t phi, Float_t dp);
     virtual ~hamcTrackOut(); 
   
     Int_t Init(Int_t ispect, hamcExpt *expt);
     Int_t Generate(hamcExpt *expt);
     Float_t GetQsq() { return qsq; };
     Float_t Getthetamin() const {return thetamin;};
     Float_t Getthetamax() const {return thetamax;};
     Float_t Getphimin() const {return phimin;};
     Float_t Getphimax() const {return phimax;};
     Int_t UpdateAtDet();

  protected:

  private: 

     void ComputePvect();
     Int_t LabToTrans();
 
// Qsq between this track and input track 'trk'
     void ComputeQsqToTrack(const hamcTrack *trk);  

     Float_t theta_central;
     Float_t thetamin, thetamax, phimin, phimax;
     Float_t tgt_mass;
     Float_t theta_iteration;
     Float_t dpp,qsq;
     Int_t which_hrs;
     Float_t xfpd,yfpd,thfpd,phfpd;
     Float_t det_dist;

     hamcTrackOut(const hamcTrackOut& phys);
     hamcTrackOut& operator=(const hamcTrackOut& phys);

     TH2F *htp;

#ifndef NODICT
ClassDef (hamcTrackOut, 0)   // Outgoing track
#endif


};

#endif



   
