//  hamcExptPREX   -- Pb Radius Experiment.
//  R. Michaels  June 2008

#include "hamcExptPREX.h"
#include "hamcExpt.h"
#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcBeam.h"
#include "hamcSpecHRS.h"
#include "hamcTgtPREX.h"
#include "hamcPhyPREX.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptPREX)
#endif

using namespace std;

hamcExptPREX::hamcExptPREX() : hamcSingles("PREX")
{
  physics = new hamcPhyPREX();
  target  = new hamcTgtPREX();
}

hamcExptPREX::~hamcExptPREX() {
}

Int_t hamcExptPREX::Init(string sfile) {

// Defaults in case the setup file is empty
  hamcSingles::SetP0(1.04);
  hamcSingles::SetTheta(5.0);

  hamcSingles::Init(sfile);

  // Defaults (can be over-ridden by input file)
  static_cast<hamcSpecHRS*>(spectrom[0])->UseCollimator();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseWarmSeptum();
  static_cast<hamcSpecHRS*>(spectrom[0])->UseMatrixTrans();

  spectrom[0]->Init(this);
  spectrom[0]->Print();

  prex_xy1 = new TH2F("prex_xy1","XY focal (main peak)",100,-1.1,0.2,100,-0.1,0.1);
  prex_xy2 = new TH2F("prex_xy2","XY focal (A_T hole)",100,-1.1,0.2,100,-0.1,0.1);
  prex_xy3 = new TH2F("prex_xy3","XY focal (A_T detector)",100,-1.1,0.2,100,-0.1,0.1);
  prex_x1  = new TH1F("prex_x1","X focal (meters)",100,-1.1,0.2);
  prex_x2  = new TH1F("prex_x2","X focal",100,-1.1,0.2);

  sumr_pc = 0;
  xcnt_pc = 0;

  return OK;
}

void hamcExptPREX::EventAnalysis() {

  hamcSingles::EventAnalysis();

  Float_t dedxX = 0.24;  // corresponding to 20 MeV, or 0.7 cm Be

  Float_t x = event->trackout[0]->xtrans;
  Float_t y = event->trackout[0]->ytrans;
  Float_t z = physics->GetCrossSection();

  if (event->inaccept == 1) {

    if (event->trackout[0]->ms_collim == 0) {

       prex_xy1->Fill(x,y,z);
       prex_x1->Fill(x,z);

    } else {

      prex_xy2->Fill(x-dedxX,y,z);
      prex_x2->Fill(x-dedxX,z);

      if (x-dedxX > -0.27 && x-dedxX < -0.23 && y > -0.02 && y < 0.003) {

        prex_xy3->Fill(x-dedxX,y,z);

        Int_t mtl_idx = target->GetMtlIndex();
        if (mtl_idx < 0 || mtl_idx >= target->GetNumMtl()) {
          cout << "hamcExptPrex::EventAna:ERROR:  bad mtl index"<<endl;
          return;
        }
 
        Float_t anum = target->GetAscatt();    // atomic num.
        if (anum == 0) {
          cout << "hamcExptPREX::EventAna:ERROR:  A = 0 ?"<<endl;
          return;
        }

         Float_t tdens = target->GetMtlDensity(mtl_idx);  // tgt density (g/cm^3)
         Float_t tlen = target->GetMtlLen(mtl_idx);  // tgt len (m)
  tlen = tlen*100;                        // need cm
         Float_t current = event->beam->beam_current;  // microAmps (uA)
         current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec

         Float_t crsec = physics->GetCrossSection();  // barns/str
         Float_t omega = 1.6e-4;   // steradians (approximate !).
         Float_t rate = current * crsec * 0.602 * tlen * tdens * omega / anum;
         sumr_pc += rate;
         xcnt_pc += 1.0;
      }
    }
  }
}

void hamcExptPREX::RunSummary() {

  cout << "hamcExptPREX:: summary "<<endl;
  cout << "Rate in A_t detector "<<
      sumr_pc/xcnt_pc<<"  Hz, for "<<xcnt_pc<<"  hits"<<endl;

  hamcSingles::RunSummary();

}
