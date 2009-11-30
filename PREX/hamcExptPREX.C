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
#include "hamcInout.h"
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
  solid_athole = 3.6e-5;  // ~1% of total solid angle.
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

  prex_xy1 = new TH2F("prex_xy1","XY at det (main peak)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy2 = new TH2F("prex_xy2","XY at det (A_T hole)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy3 = new TH2F("prex_xy3","XY at det (A_T and in Y band)",
            100,-0.7,0.2,100,-0.1,0.1);
  prex_xy4 = new TH2F("prex_xy4","XY in main detector",100,-0.7,0.2,100,-0.1,0.1);

  prex_x1  = new TH1F("prex_x1","X at det (not in A_T)",200,-0.7,0.2);
  prex_x2  = new TH1F("prex_x2","X at det (from A_T)",200,-0.7,0.2);
  prex_x3  = new TH1F("prex_x3","X at det (A_T and Y band)",
            200,-0.7,0.2);

  prex_theta = new TH1F("prex_theta","Theta generated",600,2,8);
  prex_indet = new TH1F("prex_indet","Theta in detector",600,2,8);

  sumr_pc1 = 0;
  xcnt_pc1 = 0;

  inout->AddToNtuple("inat",&inatdet);

  return OK;
}

void hamcExptPREX::EventAnalysis() {

  hamcSingles::EventAnalysis();

// Update Mar 14, 2009: there is only 1 A_T hole at the moment.  (X1)
// Update Mar 21. The dE/dX loss is done in the event processing now.

  Float_t atrl = GetSpectrom(0)->collim2_radlen1; // A_T hole rad.len.

// Cuts to define main detector

  Float_t xmain_lo, xmain_hi, ymain_lo, ymain_hi;
  xmain_lo = -0.07;
  xmain_hi = 0.04;
  ymain_lo = -0.06;
  ymain_hi = 0.07;


// Cuts to define location of A_T detector

  Float_t xmid;
  if (atrl > 0.04 && atrl < 0.045) {
    xmid = -0.095;
    xdetlo = xmid - 0.04;
    xdethi = xmid + 0.04;
  } else { 
    xmid = -0.150;
    xdetlo = xmid - 0.045;
    xdethi = xmid + 0.045;
  }

  ydetlo = -0.042;  // can be -0.042 or -0.019
  ydethi = 0.004;   // can be -0.019 or +0.004

  Float_t x = event->trackout[0]->xtrans;
  Float_t th = event->trackout[0]->thtrans;
  Float_t y = event->trackout[0]->ytrans;
  Float_t ph = event->trackout[0]->phtrans;
  Float_t z = physics->GetCrossSection();

  Float_t th_deg = (180.0/PI) * event->trackout[0]->GetTheta();

  prex_theta->Fill(th_deg);

// Extrapolate to Z of A_T detector (42.5 cm from focal plane)
  Float_t zextr = 0.425;
  Float_t xextr = x + th*zextr;
  Float_t yextr = y + ph*zextr;

  inatdet = 0;

  if (event->inaccept == 1) {

    if (event->trackout[0]->ms_collim == 0) {

       prex_xy1->Fill(xextr,yextr,z);

       if (yextr > ydetlo && yextr < ydethi) prex_x1->Fill(xextr,z);  // All within Y band

    } else {

      if( event->trackout[0]->ms_collim == 1) {
           prex_xy2->Fill(xextr,yextr,z);
           prex_x2->Fill(xextr,z);
      }


// A_T detector
      if (xextr > xdetlo && xextr < xdethi && yextr > ydetlo && yextr < ydethi) {

        inatdet = 1;

        prex_xy3->Fill(xextr,yextr,z);
        if (event->trackout[0]->ms_collim==1) prex_x3->Fill(xextr,z);  // from A_T hole and in det.

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

    }

    // In focal plane, now demand in main detector

    if (xextr > xmain_lo && xextr < xmain_hi && 
        yextr > ymain_lo && yextr < ymain_hi) {

         prex_xy4->Fill(xextr,yextr,z);
         prex_indet->Fill(th_deg);

    }



  }

}
void hamcExptPREX::RunSummary(Int_t iteration) {

  cout << "hamcExptPREX:: summary (v4.1) "<<endl;
  cout << "iteration "<<iteration<<"   "<<numiter<<endl;
  cout << "A_t detector RL "<<GetSpectrom(0)->collim2_radlen1<<endl;
  cout << "A_t detector boundary "<<xdetlo<<"  "<<xdethi<<"  "<<ydetlo<<"  "<<ydethi<<endl;
  if(xcnt_pc1 != 0) {
     cout << "Rate in A_t detector  "<<
     sumr_pc1/xcnt_pc1<<"  Hz, for "<<xcnt_pc1<<"  hits"<<endl;
  } else {
    cout << "No hits in A_t det"<<endl;
  }
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
          acc[imodel*num_mtl + idx]->Print();
      }

      sum_asy += rate * asy;
      sum_rate += rate;

// Here it's assumed diamond gets subtracted.
// Evaluate sensitivity for lead (only) 

      if ( (target->GetMtlName(idx) == "lead") ||
            (target->GetMtlName(idx) == "calcium") ||
  	     (target->GetMtlName(idx) == "tin") ) {

        if (imodel==0) asy0 = asy;
// For 2 HRS:
        xcnt = 2.0 * rate * run_time * 3600;
        asy_err = 0;   
        if (xcnt != 0) asy_err = 1e6 / TMath::Sqrt(xcnt);
        daa = asy_err / rawasy;
        cout << endl<<"+++++++++++++  MAIN TARGET +++++++++++++"<<endl;
        cout << "Material = "<<target->GetMtlName(idx)<<endl;
        cout << "Rate in 1 HRS "<<rate<<endl;
        cout << "Num of counts (sum of 2 HRS)  "<<xcnt<<endl;
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
          cout << "Sensitivity =  "<<100*sensi<<" %"<<endl;
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
 // For 2 HRS
       xcnt = 2.0 * sum_rate * run_time * 3600;
       asy_err = 0;   
       if (xcnt != 0) asy_err = 1e6 / TMath::Sqrt(xcnt);
       daa = asy_err / avg_rawasy;
       cout << endl << "Material - Averaged Results "<<endl;
       cout << "<A>_phy = "<<avg_asy;
       cout << endl << "Raw measured <A>_raw = "<<avg_rawasy;
       cout << " +/- "<<asy_err<<"   ppm "<<endl;
       cout << "Num of counts (2 HRS summed) : "<<xcnt<<endl;
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
