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
#include "hamcAccAvg.h"
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
  dpp_cut = 0.003;
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
  prex_x2  = new TH1F("prex_x2","X focal (meters)",100,-1.1,0.2);

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


  Float_t sum_rate, sum_asy;
  Float_t rate, asy, avg_asy;
  Float_t xcnt, asy_err;
  Float_t omega;
  Float_t pol = event->beam->polarization;
  Float_t polerr = event->beam->polerr;
  Float_t asy0 = 0;  // unstretched asy
  Float_t asy1 = 0;  // stretched R_N asy
  Float_t daa, drr, sensi, blowup, drrtot;

  for (Int_t imodel=0; imodel<num_phyt; imodel++) {

    cout << endl << "model  "<<imodel<<"  ----------- "<<endl;

    sum_rate = 0;
    sum_asy = 0;

    for (Int_t idx = 0; idx < num_mtl; idx++) {

      acc[imodel*num_mtl + idx]->RunSummary();

      rate = acc[imodel*num_mtl + idx]->GetRate();
      asy  = acc[imodel*num_mtl + idx]->GetAsy();
      omega = acc[imodel*num_mtl + idx]->GetOmega();

      xcnt = rate * run_time * 3600;  // run_time was in hours

      if (xcnt == 0) {
  	  cout << "hamcExptPREX::RunSum::WARNING: no hits for mtl = ";
          cout << idx<<"    model = "<<imodel<<endl;
      } else {
          cout << "\nMaterial "<<idx<<"  "<<target->GetMtlName(idx)<<endl;
          cout << "Rate "<<rate<<" Hz "<<endl;
          cout << "<A>_phys = "<<asy<<endl;
          cout << "<A>_raw = "<<asy*pol<<endl;
          cout << "omega "<<omega<<"  str "<<endl;
      }

      sum_asy += rate * asy;
      sum_rate += rate;

    }

    if (sum_rate == 0) {
      cout << "\nhamcExptPREX::RunSum::ERROR: no summed rate ?"<<endl;
    } else {
      avg_asy = sum_asy / sum_rate;  // physics asymmetry
      avg_asy = avg_asy * pol;       // raw asymmetry
      xcnt = sum_rate * run_time * 3600;
      asy_err = 0;   
      if (xcnt != 0) asy_err = 1e6 / TMath::Sqrt(xcnt);
      daa = asy_err / avg_asy;

      if (num_phyt > 1) cout << "Physics model "<<imodel<<endl;
      cout << endl << "Raw measured <A>_raw = "<<avg_asy;
      cout << " +/- "<<asy_err<<"   ppm "<<endl;
      if (imodel == 0) {
        cout << "Num of counts "<<xcnt<<endl;
        cout << "Stat precision "<<daa<<endl;
        cout << "using polarization = "<<pol<<endl;
        cout << "Total rate "<<sum_rate<<"  Hz "<<endl;
        cout << "run time "<<run_time<<" hours "<<endl;
        cout << "beam current "<<event->beam->beam_current<<" uA"<<endl; 
      }
      if (imodel==0) asy0 = avg_asy;
      if (imodel==1) {
          asy1 = avg_asy;
          sensi = 0;
          if (asy0 != 0) {
            sensi = (asy1-asy0)/asy0;          
            if (sensi < 0) sensi = -1.0*sensi;
	  }
          drr = 999999;
          if (sensi != 0) drr = 0.01 * daa / sensi;
          blowup = 1;
          if (daa != 0) blowup = (TMath::Sqrt(polerr*polerr + daa*daa))/daa;
          drrtot = drr * blowup;
          cout << endl << "Sensitivity "<<sensi<<endl;
          cout << "dR/R = "<<drr<<endl;
          cout << "blowup factor "<<blowup<<endl;
          cout << "total dR/R = "<<drrtot<<endl;
      }

    }

  }

  hamcExpt::RunSummary();


}
