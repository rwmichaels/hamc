//  hamcSingles   -- base class for single-arm experiment
//  R. Michaels  Sept 2007

#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcSpecHRS.h"
#include "hamcTarget.h"
#include "hamcBeam.h"
#include "hamcPhysics.h"
#include "hamcKine.h"
#include "hamcInout.h"
#include "THaString.h"
#include "Rtypes.h"
#include "TMath.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcSingles)
#endif


hamcSingles::hamcSingles(string sname) : hamcExpt(sname),num_mtl(0)
{
  event = new hamcEvent();  // should perhaps be a SinglesEvent ?
}

hamcSingles::~hamcSingles() {
   if (event) delete event;
}


Int_t hamcSingles::Init(string sfile) {

  hamcExpt::InitInput(sfile);
  THaString strin;
  Int_t which = RIGHTHRS; // default

  vector<string> sdata; 
  sdata = inout->GetStrVect("HRS_arm");
  if (sdata.size()>=1) {
    strin = sdata[0];
    if (strin.CmpNoCase("left")==0) {
      cout << "HRS arm = left "<<endl;
      which = LEFTHRS;
    }
    if (strin.CmpNoCase("right")==0) {
      cout << "HRS arm = right "<<endl;
      which = RIGHTHRS;
    }
  }
  sdata = inout->GetStrVect("HRS_P0");
  if (sdata.size()>=1) {
     sscanf(sdata[0].c_str(),"%f",&P0);
  }
  sdata = inout->GetStrVect("HRS_angle");
  if (sdata.size()>=1) {
     sscanf(sdata[0].c_str(),"%f",&angle);
  }

  cout << "hamcSingles:  P0 = "<<P0<<"    angle = "<<angle<<endl;

  SetSpectrom(which, P0, angle);

  hamcExpt::Init(sfile);

  num_mtl = target->GetNumMtl();

  cout << "hamcSingles::Number of target materials "<<num_mtl<<endl; 

  sumasy  = new Float_t[num_mtl];
  sumrate = new Float_t[num_mtl];
  xevtcnt = new Float_t[num_mtl];

  for (Int_t i=0; i<num_mtl; i++) {
    sumasy[i] = 0;
    sumrate[i] = 0;
    xevtcnt[i] = 0;
  }

  htpa1 = new TH2F("htpa1","Theta-Phi accepted",
              100,-1,1,100,-1,1);

  htpa2 = new TH2F("htpa2","Theta-Phi accepted",
              100,0,0.2,100,-0.62,0.62);



  return OK;

}

Int_t hamcSingles::SetSpectrom(Int_t which, Float_t pmom, Float_t theta) {

  spectrom.push_back(new hamcSpecHRS(which, pmom, theta));

  return OK;
}

void hamcSingles::EventAnalysis() {

  Int_t ldebug = 0;

  if (ldebug==2) {
    cout << "Event analysis"<<endl;
    cout << "acceptance flag "<<event->inaccept<<endl;
  }

  if ( !event->inaccept ) return;  // not in acceptance.

  Int_t mtl_idx = target->GetMtlIndex();
  if (mtl_idx < 0 || mtl_idx >= target->GetNumMtl()) {
    cout << "hamcSingles::EventAna:ERROR:  bad mtl index"<<endl;
    return;
  }
 
  Float_t anum = target->GetAscatt();    // atomic num.
  if (anum == 0) {
    cout << "hamcSingles::EventAna:ERROR:  A = 0 ?"<<endl;
    return;
  }

  Float_t domega = physics->kine->acell->domega;
  physics->kine->IncrementAcceptance();

  Float_t tdens = target->GetMtlDensity(mtl_idx);  // tgt density (g/cm^3)

  Float_t tlen = target->GetMtlLen(mtl_idx);  // tgt len (m)
  tlen = tlen*100;  // cm

  Float_t current = event->beam->beam_current;  
  current = current * 6.25e12;  // 100 uA = 6.25e14 e- / sec

  Float_t crsec = physics->GetCrossSection();  // barns/str
  Float_t asy = 1e6 * physics->GetAsymmetry();

  Float_t rate = 
       current * crsec * 0.602 * domega * tlen * tdens / anum;

  if (ldebug) {
    cout << "Singles event analysis "<<endl;
    cout << "mtl_idx "<<mtl_idx<<"  num mtl "<<num_mtl<<"   anum "<<anum<<endl;
    cout << "tgt len "<<tlen<<"   density "<<tdens<<endl;
    cout << "solid angle "<<domega<<"   beam "<<current<<endl;
    cout << "crsec "<<crsec<<"  barns/str "<<endl;
    cout << "asymmetry "<<asy<<"  ppm "<<endl;
    cout << "rate  "<<rate<<endl;
  }

  if (mtl_idx >= 0 && mtl_idx < num_mtl) {
    xevtcnt[mtl_idx] += 1;
    sumasy[mtl_idx]  += rate*asy;
    sumrate[mtl_idx] += rate;
  }

}

void hamcSingles::RunSummary() {

  Float_t sum_rate, sum_asy;
  Float_t rate, asy, avg_asy;
  Float_t xcnt, asy_err;

  sum_rate = 0;
  sum_asy = 0;

  for (Int_t idx = 0; idx < num_mtl; idx++) {

    if (xevtcnt[idx] == 0 || sumrate[idx] == 0) {

      cout << "hamcSingles::RunSum::ERROR: no counts ?"<<endl;
      cout << idx <<"  "<<xevtcnt[idx]<<"  "<<sumrate[idx]<<endl;

    } else {

      rate = sumrate[idx] / xevtcnt[idx];
      asy = sumasy[idx] / sumrate[idx];
      xcnt = rate * run_time * 3600;  // run_time in hours

      if (xcnt == 0) {
	cout << "hamcSingles::RunSum::ERROR: no run time ?"<<endl;
      } else {
        cout << "Material "<<idx<<"  "<<target->GetMtlName(idx)<<endl;
        cout << "rate "<<rate<<" Hz    num cnts "<<xcnt<<endl;
        cout << "<A> = "<<asy<<"  ppm "<<endl;
      }

      sum_asy += rate * asy;
      sum_rate += rate;

    }

  }

  if (sum_rate == 0) {
    cout << "hamcSingles::RunSum::ERROR: no summed rate ?"<<endl;
  } else {
    avg_asy = sum_asy / sum_rate;
    xcnt = sum_rate * run_time * 3600;
    asy_err = 1e6 / TMath::Sqrt(xcnt);
    cout << endl << "Overall <A> = "<<avg_asy;
    cout << " +/- "<<asy_err<<"   ppm "<<endl;
    cout << "Total rate "<<sum_rate<<"  Hz "<<endl;
    cout << "run time "<<run_time<<" hours "<<endl;
  }

// Check acceptance model

  Int_t numcell = physics->kine->acell->numcell;

  Float_t th,ph;
  Float_t domega = 0;

  for (Int_t ix = 0; ix < numcell; ix++) {
 
    for (Int_t iy = 0; iy < numcell; iy++) {

      // Put into histogram if num > cut

      Float_t cell_cut = 1;

      Float_t xcnt = physics->kine->acell->Num(ix*numcell+iy);

      if (xcnt  > cell_cut) {

        th = physics->kine->acell->GetTheta(ix);
        ph = (PI/2) - physics->kine->acell->GetPhi(iy);

	//        cout << "Xcnt "<<ix<<"  "<<iy<<"  "<<th<<"  "<<ph<<"  "<<xcnt<<endl;

	domega += physics->kine->acell->domega;

        htpa1->Fill(th,ph,xcnt);
        htpa2->Fill(th,ph,xcnt);

      }
    }
  }

  // problems: no cell cnt (Num=0), mtl_idx .ne.2

  cout << "Total solid angle "<<domega<<"  str "<<endl;

  hamcExpt::RunSummary();

}
   
