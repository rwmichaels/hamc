//  hamcExptPREX   -- Pb Radius Experiment.
//  R. Michaels  June 2008

#include "hamcExptPREX.h"
#include "hamcExpt.h"
#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcKine.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcBeam.h"
#include "hamcSpecHRS.h"
#include "hamcAccAvg.h"
#include "hamcTgtPREX.h"
#include "hamcPhyPREX.h"
#include "hamcInout.h"
#include "THaString.h"
#include "Rtypes.h"
#include <string>
#include <vector>

#ifndef NODICT
ClassImp(hamcExptPREX)
#endif

using namespace std;

 hamcExptPREX::hamcExptPREX() : hamcSingles("PREX")
{
  dpp_cut = 0.006;
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

  prex_th0 = new TH1F("prex_th0","Theta generated and in fp (unweighted)", 600,2.0,8.0);
  prex_th1 = new TH1F("prex_th1","Theta generated and in detector (unweighted)", 600,2.0,8.0);

  prex_xy0 = new TH2F("prex_xy0","XY in focal plane",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy1 = new TH2F("prex_xy1","XY at det (main peak)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy2 = new TH2F("prex_xy2","XY at det (A_T hole)",100,-0.7,0.2,100,-0.1,0.1);
  prex_xy3 = new TH2F("prex_xy3","XY at det (A_T and in Y band)",
            100,-0.7,0.2,100,-0.1,0.1);
  prex_xy4 = new TH2F("prex_xy4","XY in main detector",100,-0.7,0.2,100,-0.1,0.1);

  prex_x1  = new TH1F("prex_x1","X at det (not in A_T)",200,-0.7,0.2);
  prex_x2  = new TH1F("prex_x2","X at det (from A_T)",200,-0.7,0.2);
  prex_x3  = new TH1F("prex_x3","X at det (A_T and Y band)",
            200,-0.7,0.2);

  // Notation:  slices are 1,2,3...  hthXY if X=phi, Y = theta
  // and if X or Y = 0 it means integrate over that.

  hqsq00 = new TH1F("hqsq00","MC Qsq for all slices",200,  0.0, 0.02);
  hth00 = new TH1F("hth00","MC theta for all slices",120,-0.067,0.067);
  hph00 = new TH1F("hph00","MC phi for all slices",120,-0.04,0.04);

  hqsq10 = new TH1F("hqsq10","MC Qsq for phi slice 1 and all th",200,  0.0, 0.02);
  hqsq20 = new TH1F("hqsq20","MC Qsq for phi slice 2 and all th",200,  0.0, 0.02);
  hqsq30 = new TH1F("hqsq30","MC Qsq for phi slice 3 and all th",200,  0.0, 0.02);
  hqsq40 = new TH1F("hqsq40","MC Qsq for phi slice 4 and all th",200,  0.0, 0.02);
  hqsq50 = new TH1F("hqsq50","MC Qsq for phi slice 5 and all th",200,  0.0, 0.02);
  hqsq60 = new TH1F("hqsq60","MC Qsq for phi slice 6 and all th",200,  0.0, 0.02);
  hqsq70 = new TH1F("hqsq70","MC Qsq for phi slice 7 and all th",200,  0.0, 0.02);

  hth10 = new TH1F("hth10","MC theta for phi slice 1",120,-0.067,0.067);
  hth20 = new TH1F("hth20","MC theta for phi slice 2",120,-0.067,0.067);
  hth30 = new TH1F("hth30","MC theta for phi slice 3",120,-0.067,0.067);
  hth40 = new TH1F("hth40","MC theta for phi slice 4",120,-0.067,0.067);
  hth50 = new TH1F("hth50","MC theta for phi slice 5",120,-0.067,0.067);
  hth60 = new TH1F("hth60","MC theta for phi slice 6",120,-0.067,0.067);
  hth70 = new TH1F("hth70","MC theta for phi slice 7",120,-0.067,0.067);

  // cross check
  hth01 = new TH1F("hth01","MC theta for th slice 1",120,-0.067,0.067);
  hth02 = new TH1F("hth02","MC theta for th slice 2",120,-0.067,0.067);
  hth03 = new TH1F("hth03","MC theta for th slice 3",120,-0.067,0.067);
  hth04 = new TH1F("hth04","MC theta for th slice 4",120,-0.067,0.067);
  hth05 = new TH1F("hth05","MC theta for th slice 5",120,-0.067,0.067);
  hth06 = new TH1F("hth06","MC theta for th slice 6",120,-0.067,0.067);
  hth07 = new TH1F("hth07","MC theta for th slice 7",120,-0.067,0.067);

  hqsq01 = new TH1F("hqsq01","MC Qsq for all phi and and th slice 1",200,  0.0, 0.02);
  hqsq02 = new TH1F("hqsq02","MC Qsq for all phi and and th slice 2",200,  0.0, 0.02);
  hqsq03 = new TH1F("hqsq03","MC Qsq for all phi and and th slice 3",200,  0.0, 0.02);
  hqsq04 = new TH1F("hqsq04","MC Qsq for all phi and and th slice 4",200,  0.0, 0.02);
  hqsq05 = new TH1F("hqsq05","MC Qsq for all phi and and th slice 5",200,  0.0, 0.02);
  hqsq06 = new TH1F("hqsq06","MC Qsq for all phi and and th slice 6",200,  0.0, 0.02);
  hqsq07 = new TH1F("hqsq07","MC Qsq for all phi and and th slice 7",200,  0.0, 0.02);

  hph01 = new TH1F("hph01","MC phi for theta slice 1",120,-0.04,0.04);
  hph02 = new TH1F("hph02","MC phi for theta slice 2",120,-0.04,0.04);
  hph03 = new TH1F("hph03","MC phi for theta slice 3",120,-0.04,0.04);
  hph04 = new TH1F("hph04","MC phi for theta slice 4",120,-0.04,0.04);
  hph05 = new TH1F("hph05","MC phi for theta slice 5",120,-0.04,0.04);
  hph06 = new TH1F("hph06","MC phi for theta slice 6",120,-0.04,0.04);
  hph07 = new TH1F("hph07","MC phi for theta slice 7",120,-0.04,0.04);

  // cross check
  hph10 = new TH1F("hph10","MC phi for phi slice 1",120,-0.04,0.04);
  hph20 = new TH1F("hph20","MC phi for phi slice 2",120,-0.04,0.04);
  hph30 = new TH1F("hph30","MC phi for phi slice 3",120,-0.04,0.04);
  hph40 = new TH1F("hph40","MC phi for phi slice 4",120,-0.04,0.04);
  hph50 = new TH1F("hph50","MC phi for phi slice 5",120,-0.04,0.04);
  hph60 = new TH1F("hph60","MC phi for phi slice 6",120,-0.04,0.04);
  hph70 = new TH1F("hph70","MC phi for phi slice 7",120,-0.04,0.04);

  hqsqmid = new TH1F("hqsqmid","MC Qsq for theta,phi in middle",200,  0.0, 0.02);

  hpx1 = new TH1F("hpx1","X at detector",200,-0.8,0.2);

  prex_theta = new TH1F("prex_theta","Theta generated",600,2,8);
  prex_indet = new TH1F("prex_indet","Theta in detector",600,2,8);

  qsqf = new TH1F("qsqf","Qsq (MC) for events in detector",200,0,0.02);

  sumr_pc1 = 0;
  xcnt_pc1 = 0;

  inout->AddToNtuple("inat",&inatdet);

  return OK;
}

void hamcExptPREX::EventAnalysis() {

  hamcSingles::EventAnalysis();

// Qsq pinhole study

  Float_t qsqloc = physics->kine->qsq;

  Float_t theta = (180./3.14159)*physics->kine->theta;
  Float_t sth = TMath::Sin(physics->kine->theta);
  Float_t thwt = 1;
  if (sth > 0) thwt = 1./sth;

  Float_t th0 = event->trackout[0]->th0;
  Float_t ph0 = event->trackout[0]->ph0;

  Float_t x = event->trackout[0]->xtrans;
  Float_t th = event->trackout[0]->thtrans;
  Float_t y = event->trackout[0]->ytrans;
  Float_t ph = event->trackout[0]->phtrans;
  Float_t wt = 1e5*physics->GetCrossSection();

  Float_t th_deg = (180.0/PI) * event->trackout[0]->GetTheta();

//  Float_t thphr = TMath::Sqrt(ph0*ph0 + th0*th0);

  if (event->inaccept == 1) {

    hpx1->Fill(x, wt);

    if (x > -0.07) { // 7 cm = 6 MeV (12.4 cm/%)

      hqsq00->Fill(qsqloc, wt);
      hth00->Fill(th0, wt);
      hph00->Fill(ph0, wt);

      if (ph0 > -0.008 && ph0 < 0.010 && th0 > -0.035 && th0 < 0.035) {
	hqsqmid->Fill(qsqloc, wt);
      }


// tg_th slices here

      if (th0 > -0.06 && th0 < -0.035) {
        hqsq01->Fill(qsqloc, wt);
        hth01->Fill(th0, wt);
        hph01->Fill(ph0, wt);
      }

      if (th0 > -0.035 && th0 < -0.025) {
        hqsq02->Fill(qsqloc, wt);
        hth02->Fill(th0, wt);
        hph02->Fill(ph0, wt);
      }

      if (th0 > -0.025 && th0 < -0.010) {
        hqsq03->Fill(qsqloc, wt);
        hth03->Fill(th0, wt);
        hph03->Fill(ph0, wt);
      }

      if (th0 > -0.010 && th0 < 0.010) {
        hqsq04->Fill(qsqloc, wt);
        hth04->Fill(th0, wt);
        hph04->Fill(ph0, wt);
      }

      if (th0 > 0.01 && th0 < 0.025) {
        hqsq05->Fill(qsqloc, wt);
        hth05->Fill(th0, wt);
        hph05->Fill(ph0, wt);
      }

      if (th0 > 0.025 && th0 < 0.035) {
        hqsq06->Fill(qsqloc, wt);
        hth06->Fill(th0, wt);
        hph06->Fill(ph0, wt);
      }

      if (th0 > 0.035 && th0 < 0.06) {
        hqsq07->Fill(qsqloc, wt);
        hth07->Fill(th0, wt);
        hph07->Fill(ph0, wt);
      }

// tg_ph slices here

      if (ph0 > -0.02 && ph0 < -0.010) {
        hqsq10->Fill(qsqloc, wt);
        hth10->Fill(th0, wt);
        hph10->Fill(ph0, wt);
      }

      if (ph0 > -0.01 && ph0 < -0.008) {
        hqsq20->Fill(qsqloc, wt);
        hth20->Fill(th0, wt);
        hph20->Fill(ph0, wt);
      }

      if (ph0 > -0.008 && ph0 < -0.004) {
        hqsq30->Fill(qsqloc, wt);
        hth30->Fill(th0, wt);
        hph30->Fill(ph0, wt);
      }

      if (ph0 > -0.004 && ph0 < 0.004) {
        hqsq40->Fill(qsqloc, wt);
        hth40->Fill(th0, wt);
        hph40->Fill(ph0, wt);
      }

      if (ph0 > 0.004 && ph0 < 0.010) {
        hqsq50->Fill(qsqloc, wt);
        hth50->Fill(th0, wt);
        hph50->Fill(ph0, wt);
      }

      if (ph0 > 0.01 && ph0 < 0.015) {
        hqsq60->Fill(qsqloc, wt);
        hth60->Fill(th0, wt);
        hph60->Fill(ph0, wt);
      }

      if (ph0 > 0.015 ) {
        hqsq70->Fill(qsqloc, wt);
        hth70->Fill(th0, wt);
        hph70->Fill(ph0, wt);
      }


    }
  }
  


// Update Mar 14, 2009: there is only 1 A_T hole at the moment.  (X1)
// Update Mar 21. The dE/dX loss is done in the event processing now.

  Float_t atrl = GetSpectrom(0)->collim2_radlen1; // A_T hole rad.len.

// Cuts to define main detector

  Float_t xmain_lo, xmain_hi, ymain_lo, ymain_hi;
  xmain_lo = -0.1;
  xmain_hi = 9999;
  ymain_lo = -9999;
  ymain_hi = 9999;


// Cuts to define location of A_T detector

  Float_t xmid;
  if (atrl > 0.08 && atrl < 0.09) {
    xmid = -0.12;
    xdetlo = xmid - 0.04;
    xdethi = xmid + 0.04;
  } else { 
    xmid = -0.150;
    xdetlo = xmid - 0.045;
    xdethi = xmid + 0.045;
  }

  ydetlo = -0.042;  // can be -0.042 or -0.019
  ydethi = 0.004;   // can be -0.019 or +0.004


  prex_theta->Fill(th_deg);

// Extrapolate to z of A_T detector (42.5 cm from focal plane)
  Float_t zextr = 0.425;
  Float_t xextr = x + th*zextr;
  Float_t yextr = y + ph*zextr;

  inatdet = 0;

  if (event->inaccept == 1) {

    prex_xy0->Fill(xextr,yextr,wt);
    prex_th0->Fill(theta, thwt);

    if (event->trackout[0]->ms_collim == 0) {

       prex_xy1->Fill(xextr,yextr,wt);

       if (yextr > ydetlo && yextr < ydethi) prex_x1->Fill(xextr,wt);  // All within Y band

    } else {

      if( event->trackout[0]->ms_collim == 1) {
           prex_xy2->Fill(xextr,yextr,wt);
           prex_x2->Fill(xextr,wt);
      }


// A_T detector
      if (xextr > xdetlo && xextr < xdethi && yextr > ydetlo && yextr < ydethi) {

        inatdet = 1;

        prex_xy3->Fill(xextr,yextr,wt);
        if (event->trackout[0]->ms_collim==1) prex_x3->Fill(xextr,wt);  // from A_T hole and in det.

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

         qsqf->Fill(qsqloc,wt);
         prex_xy4->Fill(xextr,yextr,wt);
         prex_th1->Fill(theta, thwt);
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
          cout << "polerr "<<polerr<<endl;
          cout << "blowup factor "<<blowup<<endl;
          cout << "total dR/R = "<<drrtot<<endl;

          Float_t qsqbl = hqsq00->GetMean();
          Float_t angcut = GetAngCut();

// "Bottom line" printout (GeV, degrees, degrees, GeV^2  MHz, ppm, 4x %)
          cout << "Bottom line printout : E, theta, theta_cut  Qsq Rate dA/A sensi dR/R(raw) dR/R(tot)"<<endl;
          cout << "                       GeV, degrees, degrees, GeV^2  MHz, ppm, 4x %)"<<endl;

	  printf("\nBL:  %4.3f   %4.3f   %4.3f   %6.5f   %5.1f   %5.4f   %4.3f   %4.3f   %4.3f   %4.3f\n",event->beam->GetE0(),GetSpectrom(0)->GetScattAngle(),angcut,qsqbl,1e-6*rate,rawasy,100*daa,100*sensi,100*drr,100*drrtot);
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
