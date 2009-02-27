//  hamcSingles   -- base class for single-arm experiment
//  R. Michaels  Sept 2007

#include "hamcSingles.h"
#include "hamcEvent.h"
#include "hamcSpecHRS.h"
#include "hamcAccAvg.h"
#include "hamcTarget.h"
#include "hamcBeam.h"
#include "hamcTrackOut.h"
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


  hamcSingles::hamcSingles(string sname) : hamcExpt(sname),num_mtl(0),dpp_cut(-1)
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

  num_phyt = physics->NumberModels();

  cout << "hamcSingles::Number of physics models "<<num_phyt<<endl;

  Float_t thmin = event->trackout[0]->Getthetamin();
  Float_t thmax = event->trackout[0]->Getthetamax();
  Float_t phmin = event->trackout[0]->Getphimin();
  Float_t phmax = event->trackout[0]->Getphimax();

  for (Int_t imodel = 0; imodel<num_phyt; imodel++) {
    for (Int_t imtl = 0; imtl<num_mtl; imtl++) {
      Float_t xangle = angle*PI/180.;
      acc.push_back(new hamcAccAvg(xangle, thmin, thmax, phmin, phmax));
      Int_t idx=((Int_t)acc.size())-1;
      acc[idx]->InitHisto(imodel*num_mtl+imtl);
    }
  }

  return OK;

}

Int_t hamcSingles::SetSpectrom(Int_t which, Float_t pmom, Float_t theta) {

  spectrom.push_back(new hamcSpecHRS(which, pmom, theta));

  return OK;
}


Int_t hamcSingles::Run(Int_t maxevent) {

  if (!didinit) {
     cout << "hamcExpt::Run:ERROR: Need to Init() first."<<endl;
     return ERROR;
  }

  for (iteration = 0; iteration < numiter; iteration++) {

    for (Int_t ievt = 0; ievt < maxevent; ievt++ ) {

      if (ievt > 0 && ((ievt%10000)==1)) cout << "event "<<ievt<<endl;
 
      if (event) event->Process(this);

    }

    RunSummary(iteration);

  }

  return OK;
  
}


void hamcSingles::EventAnalysis() {

  Int_t ldebug = 0;

  if (ldebug==2) {
    cout << "Event analysis"<<endl;
    cout << "acceptance flag "<<event->inaccept<<endl;
  }

  if ( !event->inaccept ) return;  // not in acceptance.

  if (dpp_cut > 0) {  // if dpp_cut < 0, there is no dpp cut.
 
    Float_t dpp = event->trackout[0]->dpptrans;
    if (dpp < 0) dpp = -1.0*dpp;
    if (dpp > dpp_cut) return;   // must be in detector

  }   

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

  Float_t theta = physics->kine->theta;   // scattering angle
  Float_t phi = physics->kine->phi;       // azimuthal angle

  Float_t tdens = target->GetMtlDensity(mtl_idx);  // tgt density (g/cm^3)

  Float_t tlen = target->GetMtlLen(mtl_idx);  // tgt len (m)
  tlen = tlen*100;                        // need cm

  Float_t current = event->beam->beam_current;  // microAmps (uA)
  current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec

// Loop over physics models or parameter sets

  for (Int_t imodel = 0; imodel < num_phyt; imodel++) {

    Float_t crsec = physics->GetCrossSection(imodel);  // barns/str
    Float_t asy = 1e6 * physics->GetAsymmetry(imodel); // ppm

    Float_t rel_rate = 
        current * crsec * 0.602 * tlen * tdens / anum;

    Int_t idx = imodel*num_mtl + mtl_idx;
    if (idx >= 0 && idx < acc.size()) {
       acc[idx]->Increment(theta, phi, rel_rate, asy);
    }

    if (ldebug) {
      cout << "\n\nSingles event analysis "<<endl;
      cout << "physics model "<<imodel<<"  mtl_idx "<<mtl_idx<<endl;
      cout << "energy "<<physics->kine->energy<<"   angle "<<physics->kine->theta;
      cout << "   qsq "<<physics->kine->qsq<<endl;
      cout << "tgt len "<<tlen<<"   density "<<tdens<<"  "<<anum<<endl;
      cout << "beam "<<current<<endl;
      cout << "crsec "<<crsec<<"  barns/str "<<endl;
      cout << "physics asymmetry (not mult. by polar.)  "<<asy<<"  ppm "<<endl;
      cout << "theta "<<theta<<"  phi  "<<phi<<"   rad "<<endl;
      cout << "Relative rate  "<<rel_rate<<endl;
    }

  }

}
 
void hamcSingles::RunSummary(Int_t iteration) {

  Float_t sum_rate, sum_asy;
  Float_t rate, asy, avg_asy;
  Float_t xcnt, asy_err;
  Float_t omega;
  Float_t pol = event->beam->polarization;

  cout << "hamcSingles::RunSummary "<<endl;
 
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
	cout << "hamcSingles::RunSum:: no counts for mtl = "<<idx<<endl<<endl;
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
      cout << "hamcSingles::RunSum::ERROR: no summed rate ?"<<endl;
    } else {
      avg_asy = sum_asy / sum_rate;  // physics asymmetry
      avg_asy = avg_asy * pol;       // raw asymmetry
      xcnt = sum_rate * run_time * 3600;
      asy_err = 0;   
      if (xcnt != 0) asy_err = 1e6 / TMath::Sqrt(xcnt);
      Float_t daa = asy_err / avg_asy;
      if (num_phyt > 1) cout << "Physics model "<<imodel<<endl;
      cout << endl << "Raw measured <A>_raw = "<<avg_asy;
      cout << " +/- "<<asy_err<<"   ppm "<<endl;
      cout << "Num of counts "<<xcnt<<endl;
      cout << "Stat precision "<<daa<<endl;
      cout << "using polarization = "<<pol<<endl;
      cout << "Total rate "<<sum_rate<<"  Hz "<<endl;
      cout << "run time "<<run_time<<" hours "<<endl;
      cout << "beam current "<<event->beam->beam_current<<" uA"<<endl; 
    }

  }

  if (iteration == numiter) hamcExpt::RunSummary();

}
   
