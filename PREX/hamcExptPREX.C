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
  solid_athole = 1.5e-5;  
  inatdet = 0;
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

  prex_xy1 = new TH2F("prex_xy1","XY focal (main peak)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy2 = new TH2F("prex_xy2","XY focal (A_T hole 1)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy3 = new TH2F("prex_xy3","XY focal (A_T hole 2)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy4 = new TH2F("prex_xy4","XY focal (A_T det 1)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy5 = new TH2F("prex_xy5","XY focal (A_T det 2)",100,-0.7,0.2,100,-0.1,0.1);

  prex_x1  = new TH1F("prex_x1","X focal (meters)",200,-0.7,0.2);
  prex_x2  = new TH1F("prex_x2","X focal det 1",200,-0.7,0.2);
  prex_x2a  = new TH1F("prex_x2a","X focal det 1 in det",200,-0.7,0.2);
  prex_x3  = new TH1F("prex_x3","X focal det 2",200,-0.7,0.2);

  sumr_pc1 = 0;
  xcnt_pc1 = 0;

  sumr_pc2 = 0;
  xcnt_pc2 = 0;

  return OK;
}

void hamcExptPREX::EventAnalysis() {

  hamcSingles::EventAnalysis();

  // Update Mar 14, 2009: there is only 1 A_T hole at the moment.  (X1)

  Float_t dedxX1 = 0.055;  // corresponding to 4.4 MeV, or 1.5 cm Be
  Float_t dedxX2 = 0.073;  // corresponding to 5.9 MeV, or 2.0 cm Be

  Float_t ydetlo = -0.04;
  Float_t ydethi = 0;

  Float_t x = event->trackout[0]->xtrans;
  Float_t y = event->trackout[0]->ytrans;
  Float_t z = 1000.*physics->GetCrossSection();

  inatdet = 0;


  if (event->inaccept == 1) {

    if (event->trackout[0]->ms_collim == 0) {

       prex_xy1->Fill(x,y,z);

       if (y > ydetlo && y < ydethi) prex_x1->Fill(x,z);  // All within Y band

    } else {

      if( event->trackout[0]->ms_collim == 1) {
           prex_xy2->Fill(x-dedxX1,y,z);
           prex_x2->Fill(x-dedxX1,z);
      }
      if( event->trackout[0]->ms_collim == 2) {
           prex_xy3->Fill(x-dedxX2,y,z);
           prex_x3->Fill(x-dedxX2,z);
      }

// A_T det 1
      if (x-dedxX1 > -0.11 && x-dedxX1 < -0.04 && y > ydetlo && y < ydethi) {

        inatdet = 1;

        prex_xy4->Fill(x-dedxX1,y,z);
        if (event->trackout[0]->ms_collim==1) prex_x2a->Fill(x-dedxX1,z);  // from A_T hole and in det.

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
         Float_t omega = solid_athole;   // steradians 
         Float_t rate = current * crsec * 0.602 * tlen * tdens * omega / anum;
         sumr_pc1 += rate;
         xcnt_pc1 += 1.0;
      }

// A_T det 2  (obsolete, Mar 14, 2009)

      if (x-dedxX2 > -0.11 && x-dedxX2 < -0.06 && y > ydetlo && y < ydethi) {

        prex_xy5->Fill(x-dedxX2,y,z);

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
         Float_t omega = solid_athole;   // steradians 
         Float_t rate = current * crsec * 0.602 * tlen * tdens * omega / anum;
         sumr_pc2 += rate;
         xcnt_pc2 += 1.0;
      }

    }
  }
}

void hamcExptPREX::RunSummary(Int_t iteration) {

  cout << "hamcExptPREX:: summary "<<endl;
  cout << "iteration "<<iteration<<"   "<<numiter<<endl;
  cout << "Rate in A_t detector 1  "<<
      sumr_pc1/xcnt_pc1<<"  Hz, for "<<xcnt_pc1<<"  hits"<<endl;
  cout << "Rate in A_t detector 2  "<<
      sumr_pc1/xcnt_pc2<<"  Hz, for "<<xcnt_pc2<<"  hits"<<endl;
  cout << "Solid angle of A_t hole "<<solid_athole<<" str "<<endl;


  Float_t sum_rate, sum_asy;
  Float_t rate, asy, avg_asy, rawasy, avg_rawasy;
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
      rawasy = pol * asy;
      omega = acc[imodel*num_mtl + idx]->GetOmega();

      xcnt = rate * run_time * 3600;  // run_time was in hours

      if (xcnt == 0) {
  	  cout << "hamcExptPREX::RunSum::WARNING: no hits for mtl = ";
          cout << idx<<"    model = "<<imodel<<endl;
      } else {
          cout << "\nMaterial "<<idx<<"  "<<target->GetMtlName(idx)<<endl;
          cout << "Rate "<<rate<<" Hz "<<endl;
          cout << "<A>_phys = "<<asy<<endl;
          cout << "<A>_raw = "<<rawasy<<endl;
          cout << "omega "<<omega<<"  str "<<endl;
      }

      sum_asy += rate * asy;
      sum_rate += rate;

// Here it's assumed diamond gets subtracted.
// Evaluate sensitivity for lead (only) 

      if (target->GetMtlName(idx) == "lead") {
        if (imodel==0) asy0 = asy;
        xcnt = rate * run_time * 3600;
        asy_err = 0;   
        if (xcnt != 0) asy_err = 1e6 / TMath::Sqrt(xcnt);
        daa = asy_err / rawasy;
        cout << endl<<"+++++++++++++  LEAD +++++++++++++"<<endl;
        cout << "Num of counts "<<xcnt<<endl;
        cout << "Stat error = "<<100*daa<<" % "<<endl;
        if (imodel==1) {
          asy1 = asy;
          sensi = 0;
          if (asy1 != 0) {
            sensi = (asy1-asy0)/asy1;          
            if (sensi < 0) sensi = -1.0*sensi;
	  }
          drr = 999999;
          if (sensi != 0) drr = 0.01 * daa / sensi;
          blowup = 1;
          if (daa != 0) blowup = (TMath::Sqrt(polerr*polerr + daa*daa))/daa;
          drrtot = drr * blowup;
          cout << "Sensitivity for lead  "<<100*sensi<<" %"<<endl;
          cout << "dR/R = "<<drr<<endl;
          cout << "blowup factor "<<blowup<<endl;
          cout << "total dR/R = "<<drrtot<<endl;
// "Bottom line" printout (GeV, degrees, GHz, ppm, 4* %)
	  printf("\nBL:  %4.3f  %2.0f  %5.4f  %5.4f  %4.3f  %4.3f  %4.3f  %4.3f\n",event->beam->GetE0(),GetSpectrom(0)->GetScattAngle(),1e-9*rate,rawasy,100*daa,100*sensi,100*drr,100*drrtot);
          cout << "++++++++++++++++++++++++++++++++++++++++++"<<endl;
	}
      }

    } // Loop over materials

    if (imodel==0) {
     if (sum_rate == 0) {
       cout << "\nhamcExptPREX::RunSum: no rate (probably low stats)"<<endl;
     } else {
       avg_asy = sum_asy / sum_rate;     // target-averaged physics asymmetry
       avg_rawasy = avg_asy * pol;       // raw asymmetry
       xcnt = sum_rate * run_time * 3600;
       asy_err = 0;   
       if (xcnt != 0) asy_err = 1e6 / TMath::Sqrt(xcnt);
       daa = asy_err / avg_rawasy;
       cout << endl << "Material - Averaged Results "<<endl;
       cout << "<A>_phy = "<<avg_asy;
       cout << endl << "Raw measured <A>_raw = "<<avg_rawasy;
       cout << " +/- "<<asy_err<<"   ppm "<<endl;
      
       cout << "Num of counts "<<xcnt<<endl;
       cout << "Stat precision "<<daa<<endl;
       cout << "using polarization = "<<pol<<endl;
       cout << "Total rate "<<sum_rate<<"  Hz "<<endl;
       cout << "run time "<<run_time<<" hours "<<endl;
       cout << "beam current "<<event->beam->beam_current<<" uA"<<endl; 
       cout << "dpp cut "<<dpp_cut<<endl;
     }
    }

  }   // Loop over models

  if (iteration+1 == numiter) hamcExpt::RunSummary(iteration);


}
