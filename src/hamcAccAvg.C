//  hamcAccAvg   -- Class to organize the averaging over acceptance.

//  One danger:  you must run hamc with enough statistics to fill the 
//  cells of the acceptance, see the cut ncellcut.
//  You can then see that the total solid angle agrees approx. with
//  the rectangular approximation.  Another way to see it is that
//  the rate "plateaus" with number of events.  Rule of thumb: 10M events
//  to get a reaonable rate estimate.  (or reduce ncellcut, or ...)

//  R. Michaels  Jan 2009


#include "hamcAccAvg.h"
#include "hamcExpt.h"
#include "hamcPhysics.h"
#include "hamcBeam.h"
#include "hamcEloss.h"
#include "hamcEvent.h"
#include "hamcTarget.h"
#include "hamcInout.h"
#include "TRandom.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifndef NODICT
ClassImp(hamcAccAvg)
#endif

hamcAccAvg::hamcAccAvg(Float_t thc, Float_t th1, Float_t th2, Float_t ph1, Float_t ph2): 
  central_angle(thc),thmin(th1),thmax(th2),phmin(ph1),phmax(ph2),did_init(kFALSE),did_summary(kFALSE),hacc(0) {
  Init();
}

hamcAccAvg::~hamcAccAvg() {
  delete [] sumcnt;
  delete [] sumrate;
  delete [] sumasy;
}

void hamcAccAvg::Init() {
    dtheta = (thmax-thmin)/((Float_t)(numcell));
    dphi = (phmax-phmin)/((Float_t)(numcell));  
    sumcnt = new Float_t[numcell*numcell];
    sumrate = new Float_t[numcell*numcell];
    sumasy = new Float_t[numcell*numcell];
    for (Int_t i=0; i<numcell*numcell; i++) { 
        sumcnt[i] = 0; 
        sumrate[i] = 0; 
        sumasy[i] = 0;
    }
    debug = 0;
    rate = 0;
    asy = 0;
    omega = 0;
    did_init = kTRUE;
    did_summary = kFALSE;
}

void hamcAccAvg::InitHisto(Int_t id) {
  char cname[80],ctitle[80];
  sprintf(cname,"acc%d",id);
  sprintf(ctitle,"Acceptance Check %d",id);
  hacc = new TH2F(cname,ctitle,100,-0.1,0.1,100,-0.1,0.1);
}

void hamcAccAvg::Increment(Float_t th, Float_t ph, Float_t relrate, Float_t asy) {
    if ( !did_init) {
      cout << "hamcAccAvg::ERROR: object not initialized !"<<endl;
      return;
    }
    Float_t domega = TMath::Sin(th) * dtheta * dphi;
    Float_t rate = relrate * domega;
    Int_t icell, jcell, ncell;
    icell = ((Int_t)((th-thmin)/dtheta));
    jcell = ((Int_t)((ph-phmin)/dphi));
    ncell = icell*numcell + jcell;
    if (ncell >= 0 || ncell < numcell*numcell) {
        sumcnt[ncell] += 1.0;
        sumrate[ncell] += rate;
        sumasy[ncell] += asy; 
    }
}

void hamcAccAvg::Print() {
    cout<<"\nhamcAccCell parameters"<<endl;
    if ( !did_init ) cout<<" object NOT initialized "<<endl;
    cout<<"phi limits "<<phmin<<" "<<phmax<<endl;
    cout<<"theta limits "<<thmin<<" "<<thmax<<endl;
    cout<<"dtheta "<<dtheta<<"  dphi "<<dphi<<endl;
    cout<<"num cells "<<numcell<<endl;
}

Float_t hamcAccAvg::GetTheta(Int_t icell) {
    if (debug) cout << "chk gettheta "<<icell<<"  "<<thmin + dtheta*((Float_t)(icell))<<endl;
    return thmin + dtheta*((Float_t)(icell));
}

Float_t hamcAccAvg::GetPhi(Int_t icell) {
    if (debug) cout << "chk getphi "<<icell<<"  "<<phmin + dphi*((Float_t)(icell))<<endl;
    return phmin + dphi*((Float_t)(icell));
}

Float_t hamcAccAvg::GetNum(Int_t icell) { 
    if (icell < 0 || icell >= numcell*numcell) return 0;
    if (debug) cout << "chk Num "<<icell<<"  "<<numcell*numcell<<"  "<<sumcnt[icell]<<endl;
    return sumcnt[icell];
}

Float_t hamcAccAvg::GetAsy() {
  if ( Check() == ERROR ) return ERROR;
  return asy;
}

Float_t hamcAccAvg::GetRate() {
  if ( Check() == ERROR ) return ERROR;
  return rate;
}

Float_t hamcAccAvg::GetOmega() {
  if ( Check() == ERROR ) return ERROR;
  return omega;
}


void hamcAccAvg::RunSummary() {

  Float_t asysum = 0;
  Float_t ratecell;
  Int_t icell;
  Int_t lprint = 1;
  rate = 0;
  asy = 0;
  omega = 0;

  if (lprint) cout << endl<<endl<<"Cell scan ---- "<<endl;

  for (Int_t ix = 0; ix < numcell; ix++) {  // theta

    for (Int_t iy = 0; iy < numcell; iy++) {  // phi

      icell = ix*numcell + iy;

      Float_t theta = GetTheta(ix);
      Float_t phi = GetPhi(iy);
      if (lprint && theta > 0.07 && theta < 0.14 &&
        phi > 1.39 && phi<1.75 ) {
        cout << "icell " <<icell<<"   theta "<<theta<<"   phi "<<phi;
        cout << "   sum "<<sumcnt[icell]<<"  "<<sumrate[icell]<<endl;
      }

      if (GetNum(icell) < ncellcut) continue;
      if (sumrate[icell] <= 0) continue;

      ratecell = sumrate[icell]/sumcnt[icell];
 
      rate += ratecell;

// Take avg asy in cell, then weight by rate
      asysum += ratecell * sumasy[icell]/sumcnt[icell]; 

      omega += TMath::Sin(GetTheta(ix)) * dtheta * dphi;

      if (hacc) {
      
        Float_t stc, ctc, ttc, sts, cts, sps, cps, uparam;
        stc = TMath::Sin(central_angle);
        ctc = TMath::Cos(central_angle);
        ttc = TMath::Tan(central_angle);
        sts = TMath::Sin(GetTheta(ix));
        cts = TMath::Cos(GetTheta(ix));
        sps = TMath::Sin(GetPhi(iy));
        cps = TMath::Cos(GetPhi(iy));
        uparam = stc*sps*sts + ctc*cts;
        if (uparam != 0) {
           Float_t tantheta_t = cps*sts / uparam;
           Float_t tp_prime = sps*sts / uparam;
           Float_t tanphi_t = (tp_prime - ttc) / (1 + tp_prime*ttc);
           hacc->Fill(tanphi_t, tantheta_t);
	}

      }

    }
  }
 
  if (rate != 0) asy = asysum / rate;

  did_summary = kTRUE;    

}
  
Int_t hamcAccAvg::Check() {
  if (!did_init) {
    cout << "hamcAccAvg::ERROR: class not initialized "<<endl;
    return ERROR;
  }
  if (!did_summary) {
    cout << "hamcAccAvg::ERROR: class summary not done "<<endl;
    return ERROR;
  }
  return OK;
}
