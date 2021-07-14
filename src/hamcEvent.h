#ifndef ROOT_hamcEvent
#define ROOT_hamcEvent

//  hamcEvent   -- class for an event
//  R. Michaels  Sept 2007

#include "Rtypes.h"
#include <vector>

class hamcBeam;
class hamcTrackOut;
class hamcExpt;
class HrsTrkCorr;

class hamcEvent {

  public:

     hamcEvent();
     virtual ~hamcEvent();

     Int_t Init(hamcExpt *ex);
     Int_t InitBeam(hamcExpt *ex);
     Int_t Process(hamcExpt* ex);
     hamcBeam* beam;
     std::vector<hamcTrackOut*> trackout;
     hamcTrackOut* tracknoms;
     hamcTrackOut* tracksieve;
     hamcTrackOut* tracksievenoms;
     HrsTrkCorr* trkcorr;
     Int_t inaccept;  // In acceptance (1) or not (0)
     Int_t brkpoint;  // which break point we're on.

private: 

     Bool_t did_init;
     static const Int_t debug=0;
     Float_t xtgt, ytgt, thtgt, phtgt;
     Float_t xcol, ycol, thcol, phcol, xsep, ysep;
     Float_t thnoms, phnoms, thsiv, phsiv, thsinoms, phsinoms;
     Int_t evnum;
 
     hamcEvent(const hamcEvent& evt);
     hamcEvent& operator=(const hamcEvent& evt);

#ifndef NODICT
ClassDef (hamcEvent, 0)   // An event
#endif

};

#endif



   
