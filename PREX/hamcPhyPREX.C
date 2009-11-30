//  hamcPhyPREX   -- class for the physics of PREX
//  D. Jaunzeikare, R. Michaels  May 2008


#include "hamcPhyPREX.h"
#include "hamcExpt.h"
#include "hamcTarget.h"
#include "hamcEvent.h"
#include "hamcBeam.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcInout.h"
#include "hamcKine.h"
#include "THaString.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h" 
#include "TROOT.h" 
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#ifndef NODICT
ClassImp(hamcPhyPREX)
#endif


hamcPhyPREX::hamcPhyPREX() : hamcPhysics()
{
  phy_name = "PREX physics";
  scatt_process = "elastic";
  whichmodel = HORPB;  // default
  do_radiate = kTRUE;
  num_models = 2;  // if =2 we're considering stretched vs unstretched R_N
}


hamcPhyPREX::~hamcPhyPREX() { }

Int_t hamcPhyPREX::SetModel(Int_t modeln) {
  if (modeln >= HORPB && modeln <= NL3M05) {
    whichmodel = modeln;
  } else {
    cout << "hamcPhyPREX::Warning: model # outside range"<<endl;
    return ERROR;
  }
  num_models = 1;
  if (modeln == HORPB) num_models = 2; // considering stretched vs unstretched R_N
  return OK;
}


Int_t hamcPhyPREX::Init(hamcExpt* expt) {

  hamcPhysics::Init(expt);

  expt->inout->BookHisto(kTRUE, kTRUE, IFOCAL, "Asy0",
	      "Asymmetry (accepted, weighted)",
                          &asy0,200,-2e-7,2e-6);

  THaString strin;
  vector<string> sdata; 
  sdata = expt->inout->GetStrVect("PREX_model");
  if (sdata.size()>=1) {
    strin = sdata[0];
    if (strin.CmpNoCase("horpb")==0) {
      cout << "Using HORPB (original) model from Horowitz "<<endl;
      cout << "With stretching of R_N by 1%"<<endl;
      whichmodel = HORPB;
      num_models = 2;
    }
    if (strin.CmpNoCase("horca")==0) {
      cout << "Using Calcium 48 model from Horowitz "<<endl;
      cout << "With stretching of R_N by 1%"<<endl;
      whichmodel = HORCA;
      num_models = 2;
    }
    if (strin.CmpNoCase("horsn")==0) {
      cout << "Using Tin 120 model from Horowitz "<<endl;
      cout << "With stretching of R_N by 1%"<<endl;
      whichmodel = HORSN;
      num_models = 2;
    }
    if (strin.CmpNoCase("si")==0) {
      cout << "Using SI model from Horowitz "<<endl;
      whichmodel = SI;
      num_models = 1;
    }
    if (strin.CmpNoCase("siii")==0) {
      cout << "Using SIII model from Horowitz "<<endl;
      whichmodel = SIII;
      num_models = 1;
    }
    if (strin.CmpNoCase("sly4")==0) {
      cout << "Using SLY4 model from Horowitz "<<endl;
      whichmodel = SLY4;
      num_models = 1;
    }
    if (strin.CmpNoCase("fsu")==0) {
      cout << "Using FSU-Gold model from Horowitz "<<endl;
      whichmodel = FSU;
      num_models = 1;
    }
    if (strin.CmpNoCase("nl3")==0) {
      cout << "Using NL3 model from Horowitz "<<endl;
      whichmodel = NL3;
      num_models = 1;
    }
    if (strin.CmpNoCase("nl3m05")==0) {
      cout << "Using NL3M05 model from Horowitz "<<endl;
      whichmodel = NL3M05;
      num_models = 1;
    }
    if (strin.CmpNoCase("nl3p06")==0) {
      cout << "Using NL3P06 model from Horowitz "<<endl;
      whichmodel = NL3P06;
      num_models = 1;
    }

  }


  LoadFiles();  // Load the lookup files
  didinit = kTRUE;


  if (histo_test) {

    hpph1 = new TH1F("hpph1","Horowitz Lookup Crsec  vs  angle",2000,3,9);
    hpph2 = new TH1F("hpph2","Mott * FF (1.05 GeV)",2000,3,9);
    hpph3 = new TH1F("hpph3","Cross Section (0.2482 GeV)",2000,17,95);
    hpph4 = new TH1F("hpph4","Cross Section (0.502 GeV)",2000,12,42);
    hpph5 = new TH2F("hpph5","Percent diff vs angle",100,3,9,100,-10,2);
    hpph6 = new TH1F("hpph6","Asymmmetry  vs  angle",2000,3,9); 

    Float_t energy[4];
    energy[0] = 1.05;
    energy[1] = 1.05;
    energy[2] = 0.2482;
    energy[3] = 0.502;

    Float_t theta_rad, theta_degr;

    for (Int_t iene = 0; iene<4; iene++) {

    for (Int_t iang=0; iang<2000; iang++) {

      theta_degr = 3 + 125*(0.5+(Float_t)iang)/2000;
      theta_rad = 3.1415926*theta_degr/180.0;

      CrossSection(energy[iene],theta_rad,0);
      Asymmetry(energy[iene],theta_rad,0);

      Float_t frecoil = 1 + (energy[iene]/195.)*(1-TMath::Cos(theta_rad));
      Float_t eprime = energy[iene]/frecoil;
      qsq = 2*energy[iene]*eprime*(1-TMath::Cos(theta_rad));
      Float_t crsec2 = CalculateCrossSection(0, energy[iene], theta_degr);

      if (iene==0) hpph1->Fill(theta_degr,crsec);
      if (iene==1) hpph2->Fill(theta_degr,crsec2);
      if (iene==2) hpph3->Fill(theta_degr,crsec2);
      if (iene==3) hpph4->Fill(theta_degr,crsec2);
      if (iene==0) hpph6->Fill(theta_degr,asy0); 
      if (iene==1 && theta_degr>3 && theta_degr<9) {
	Float_t diff = crsec2 - crsec;
        if (crsec != 0) {
	  diff = 100*diff/crsec;
	} else {
          diff = -5;
	}
	//        cout << "crsec "<<"  "<<theta_degr<<"  "<<crsec<<"  "<<crsec2<<"  "<<diff<<endl;

        hpph5->Fill(theta_degr,diff);
      }

    }
    }

  }

  if (accept_test) {
    FILE *fd;
    char schar[200]; 
    fd=fopen("./PREX/accept.dat", "r");
    Double_t weisum=0;
    Double_t crcsum=0;
    Double_t asyavg;    
    if (fd==NULL) {
      printf("hamcPhyPREX:: trying to test acceptance model.\n");
      printf("  but cant find input.\n");
    } else {
      Float_t theta_degr, ene, xacc;
      ene = 1.05;  // GeV 
      while(fgets(schar,100,fd)!=NULL) {
        sscanf(schar, "%f  %f", &theta_degr, &xacc);

        theta_rad = 3.1415926*theta_degr/180.0;
        CrossSection(ene,theta_rad,0);
        Asymmetry(ene,theta_rad,0);

        weisum += xacc * crsec * asy0 * 1e6;
        crcsum += xacc * crsec;        

	//        cout << "Check accept "<<theta_degr<<"  "<<xacc<<"  "<<crsec<<"  "<<asy0<<endl;
       

      }

      asyavg = 0; 
      if (crcsum > 0) asyavg = weisum / crcsum;
      printf("Acc-avg asy (approx) = %7.5f \n",asyavg);

    }

  }

  if (quick_check) { // code exists after this

    Float_t theta_degr, theta_rad, ene;

    for (Int_t itry=0; itry<10; itry++) {

      theta_degr = 5.0;
      ene = 1.05;

      if (itry==1) {
	theta_degr = 3.51;
        ene = 1.2;
      }
      if (itry==2) {
	theta_degr = 4.0;
        ene = 1.2;
      }
      if (itry==3) {
	theta_degr = 4.6;
        ene = 1.25;
      }
      if (itry==4) {
	theta_degr = 4.2;
        ene = 1.3;
      }
      if (itry==5) {
	theta_degr = 4.6;
        ene = 1.3;
      }
      if (itry==6) {
	theta_degr = 5;
        ene = 1.4;
      }
      if (itry==7) {
	theta_degr = 3.8;
        ene = 1.2;
      }
      if (itry==8) {
	theta_degr = 3.8;
        ene = 1.1;
      }
      if (itry==9) {
	theta_degr = 4.2;
        ene = 1.0;
      }

      theta_rad = 3.1415926*theta_degr/180.0;
      CrossSection(ene,theta_rad,0);
      Asymmetry(ene,theta_rad,0);
      Asymmetry(ene,theta_rad,1);

      cout << "CHECK :  energy "<<ene<<"  theta "<<theta_degr<<endl;
      cout << "crsec  "<<crsec<<"    A0= "<<asy0<<"  A1= "<<asy1<<endl;
      cout << endl << "-------------------------------------------- "<<endl;

    }

    exit(0);
  }



  if (quick_feasibility) {  // code exits after this

    Int_t itgt = 0;
    if (whichmodel == HORPB) itgt = 1;

    Float_t tdens = expt->target->GetMtlDensity(itgt);  // tgt density (g/cm^3)
    Float_t tlen = expt->target->GetMtlEffLen(itgt);  // tgt len (m)
    tlen = tlen*100;                        // need cm
    Float_t current = expt->event->beam->beam_current;  // microAmps (uA)
    current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec
    Float_t anum = expt->target->GetAscatt();    // atomic num.

    cout << "Check feasibility "<<endl;
    cout << "Tgt  A = "<<anum<<"    tdens = "<<tdens<<"   tlen = "<<tlen<<endl;
    cout << "Current = "<<current<<endl;

    FILE *fd;
    char schar[200]; 
    fd=fopen("./PREX/accept.dat", "r");
    Double_t weisum=0;
    Double_t crcsum=0;
    Double_t asyavg;    
    if (fd==NULL) {
      printf("hamcPhyPREX:: Acceptance function not found !\n");
      exit(0);
    } else {

      Float_t theta_degr, ene, xacc;
      Float_t daa, rate, asyavg, domega, fom;
      Float_t daasum, ratesum, asysum, xnorm;
      domega = 0.0037;

      for (Int_t iene=0; iene < 30; iene++) {
 
        ene = 0.6 + 0.05*((Float_t)iene); // GeV

        ratesum = 0;
        xnorm=0;
        asysum = 0;
        daasum = 0;

        fclose(fd);
        fd=fopen("./PREX/accept.dat", "r");

        while(fgets(schar,100,fd)!=NULL) {
          sscanf(schar, "%f  %f", &theta_degr, &xacc);

          theta_rad = 3.1415926*theta_degr/180.0;
          CrossSection(ene,theta_rad,0);
          Asymmetry(ene,theta_rad,0);
          Asymmetry(ene,theta_rad,1);
          daa = 0;
   	  if (asy0 != 0) daa = (asy1-asy0)/asy0;
          rate = xacc * current * crsec * 0.602 * tlen * tdens * domega / anum;
          ratesum += rate;  
          daasum += daa * rate;
          asysum += asy0 * rate;
          xnorm += xacc;
	  if (ene > 1.2 && ene < 1.1) { 
	     cout << "energy "<<ene<<"  theta "<<theta_degr<<endl;
	     cout << "accept  "<<xacc<<"  A= "<<asy0<<"  "<<asy1<<"  daa = "<<daa<<endl;
	     cout << "crsec  "<<crsec<<endl;
	     cout << "rate "<<rate<<"  "<<ratesum<<"  cnt = "<<xnorm<<endl;
	  }
	}
        if (ratesum > 0) {        
          asyavg = asysum / ratesum;
          daa = daasum / ratesum;
          rate = ratesum / xnorm;
          fom = 1e9 * rate * (asyavg*asyavg) * (daa*daa);
	} else {
          asyavg=0; daa=0; rate=0; fom=0;
	}
        cout << "E= "<<ene<<"   <A> = "<<1e6*asyavg<<" ppm";
        cout << "  rate = "<<rate<<" Hz     daa = "<<daa<<"   fom = "<<fom<<endl;
        

      }
    }

    exit(0);   // done

  }  // quick_feasibility
    

  return 1;
}


Int_t hamcPhyPREX::Generate(hamcExpt *expt) {

   Int_t ldebug=0;
   static int lentry=0;

   Float_t energy = kine->energy;    // GeV
   Float_t theta = kine->theta;      // radians
   qsq = kine->qsq;                  // GeV^2

   Float_t anum = expt->target->GetAscatt();

   Int_t mtl_idx = expt->target->GetMtlIndex();
   Float_t tdens = expt->target->GetMtlDensity(mtl_idx);  // tgt density (g/cm^3)
   Float_t tlen = expt->target->GetMtlLen(mtl_idx);  // tgt len (m)
   Float_t tefflen = expt->target->GetMtlEffLen(mtl_idx);  // eff. tgt len (m)

   if (anum == 12) {
     crsec = CalculateCrossSection(1, energy, theta*180/PI);
     asymmetry = CalculateAsymmetry(1);
     drate = CalculateDrate(anum, tdens, tlen, crsec);
     if (ldebug) {
       cout << endl<< "C12:  energy "<<energy<<"  theta "<<theta<<"  qsq "<<qsq<<endl;
       cout<< "C12 crsec = "<<crsec<<"   A = "<<asymmetry<<endl<<endl;
     }
     return OK;
   }

// Compute the cross section for PREX
   CrossSection(energy, theta, 0);

// Compute the PV asymmetry 
   if (num_models > 1) Asymmetry(energy, theta,1);  // stretched
   Asymmetry(energy, theta,0);  // unstretched (call this last)

   Drate(anum, tdens, tefflen, crsec); //Hz/uA

   return OK;
}

Float_t hamcPhyPREX::GetCrossSection(Int_t istretch) const {
// No model dependence yet (Jan 2009)

  return crsec;

}


Float_t hamcPhyPREX::GetAsymmetry(Int_t istretch) const {

  if (istretch == 1) return asy1;
  return asy0;

}



Int_t hamcPhyPREX::CrossSection(Float_t energy, Float_t angle_rad, Int_t stretch) {

  Float_t angle = angle_rad * 180 / PI;

  crsec = 0;   

  if (!didinit) {
    cout << "ERROR: hamcPhyPREX:: Not initialized."<<endl;
    return -1;
  }

  Int_t debug = 0;
  // Lookup the cross section for this E,theta(angle)

  /*find the index of angle and energy in the class variables angle_row and energy_row respectively */
  Int_t indxAngle = FindAngleIndex(angle); 

  if (indxAngle <= 0) return -1; 

  Int_t indxEnergy = FindEnergyIndex(energy);

  if (indxEnergy <= 0) return -1;

  Float_t crsc1, crsc2;

  //find the angles above and below the actual angle; used for interpolation
  Float_t angle_upper =angle_row[indxAngle];
  Float_t angle_lower = angle_row[indxAngle-1];

  crsc1 = crsc_tables[stretch][indxEnergy][indxAngle]; /*get the cross section value for angle value larger than the actual*/
  crsc2 = crsc_tables[stretch][indxEnergy][indxAngle+1]; //get the cross section value for energy value one below the actual
    // }

  crsec = Interpolate(angle_lower, angle_upper, angle, crsc1, crsc2)/1000;  

  if (debug) {
     cout << "\n LEAD Cross section : "<<endl;
     cout << "stretch "<<stretch<<"   indices "<<indxEnergy<<"   "<<indxAngle+1<<endl;
     cout << "energy " << energy << " GeV " << "angle " << angle << endl;
     cout <<"Interpolation of crsec "<<crsc1/1000 <<" "<<crsc2/1000<<" "<<crsec<<endl;
  }

  Float_t calccrsec = CalculateCrossSection(0, energy, angle);
  if (debug) cout <<"calculated lead cross section " << calccrsec <<endl;

  return OK;
}


Int_t hamcPhyPREX::Asymmetry(Float_t energy, Float_t angle_rad, Int_t stretch) {

  if (!didinit) {
    cout << "ERROR: hamcPhyPREX:: Not initialized."<<endl;
       return -1;
  }
  // lookup the asymmetry for this E,theta

  Float_t angle = angle_rad * 180 / PI;

  asymmetry = -9999;

  Int_t indxAngle = FindAngleIndex(angle);
  Int_t indxEnergy = FindEnergyIndex(energy);

  if (indxAngle <= 0 || indxEnergy <= 0) return -1;

  Float_t asymmetry1, asymmetry2;
  Float_t angle_upper =angle_row[indxAngle] ;
  Float_t angle_lower = angle_row[indxAngle-1];

  asymmetry1 = asymmetry_tables[stretch][indxEnergy][indxAngle];
  asymmetry2 = asymmetry_tables[stretch][indxEnergy][indxAngle+1];
  asymmetry = Interpolate(angle_lower, angle_upper, angle, asymmetry1, asymmetry2);

  if (stretch == 0) asy0 = asymmetry;
  if (stretch == 1) asy1 = asymmetry;

  int debug=0;

  if (debug) {
    cout<<"angle  "<<angle<<endl;
    cout << "Indices "<<indxAngle<<"  "<<indxEnergy<<"   "<<stretch<<endl;
    cout << "Lead Asymmetries  "<<asymmetry1 << " "<<asymmetry2 <<" " <<asymmetry<<endl;   
    //    Float_t calcasy = CalculateAsymmetry(0);
    //    cout << "Calculated lead asy "<<calcasy<<endl;
    cout << "stretch "<<stretch<<"  "<<asy0<<"  "<<asy1<<endl;
    cout << "size = "<<asymmetry_tables.size()<<endl;
  }

  return 1;
}

Int_t hamcPhyPREX::Drate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec) {

  tlen = tlen*100;   // need cm
  
  //cout<<"tlen="<<tlen<<", tdens = "<<tdens<<", anum="<<anum<<", crsec = "<<crsec<<endl;
  Float_t avg_omega = 0.004671;
  drate = 6.25e12 * crsec * 0.602 * tlen * tdens * avg_omega / anum; //Hz/uAww
  return 1;
}


Int_t hamcPhyPREX::FindAngleIndex(Float_t angle){

  vector<Float_t>::const_iterator iterAngle = upper_bound(angle_row.begin(), angle_row.end(), angle); //search for the first value of angle which is larger than actual

  return iterAngle - angle_row.begin();
  
}

Int_t hamcPhyPREX::FindEnergyIndex(Float_t energy){
  
  // Adjust energy up 4.5% (the acceptance) to be closer to right value
  // and to not over-estimate the rate.  (Problem is 50 MeV steps)
  Float_t rounded_energy = round(1.045*energy*20)/20;
  
  vector<Float_t>::const_iterator iterEnergy  = lower_bound(energy_row.begin(), energy_row.end(), rounded_energy); //search for the first value of energy which is larger or equal than actual 

  UInt_t isize = energy_row.size();
  
  if (rounded_energy > energy_row[isize-1]) return -1;

  return iterEnergy - energy_row.begin()-1; //indexing in vectors starts from 0, therefore -1
}

Float_t hamcPhyPREX::Interpolate(Float_t min, Float_t max, Float_t mid, Float_t val1, Float_t val2){
   /*This function interpolates linearly between two values (val1 and val2)
  for example: Interpolate(angle_lower, angle_upper, angle, crsc1, crsc2)
  */
  
  Float_t RANGE = max - min;

  Float_t per1 = (max - mid) / RANGE;
  Float_t per2 = (mid - min) / RANGE;

  Float_t val = val1*per1 + val2*per2;

  //  int debug=0;
  
  return val;
}

Int_t hamcPhyPREX::LoadFiles() {
  /* Loading Horowitz tables in memory. Data is saved in 3D vector( aVector[stretch][energy][angle]*/

  int stretch = 0;

  vector<vector<Float_t> > crsc_table_temp;
  vector<vector<Float_t> > asymmetry_table_temp;

  LoadHorowitzTable(crsc_table_temp, asymmetry_table_temp,stretch); 
  crsc_tables.push_back(crsc_table_temp);
  asymmetry_tables.push_back(asymmetry_table_temp);

  crsc_table_temp.clear();
  asymmetry_table_temp.clear();

  if (whichmodel == HORPB ||
      whichmodel == HORCA ||
      whichmodel == HORSN ) {  
    stretch = 1; 
    LoadHorowitzTable(crsc_table_temp, asymmetry_table_temp,stretch);
    crsc_tables.push_back(crsc_table_temp);
    asymmetry_tables.push_back(asymmetry_table_temp);
  }

  LoadFormFactorTable();

  LoadC12FormFactorTable();
  
  cout<< "hamcPhyPREX:  Tables loaded" <<endl;


  return 1;
}


Float_t hamcPhyPREX::CalculateCrossSection(Int_t nuc, Float_t energy, Float_t angle) {

  // nuc = 0  --> lead
  // nuc = 1  --> C12
  // no other choices !

  // energy in GeV,  angle in degrees.

  if (nuc != 0 && nuc != 1) {
    cout << "hamcPhyPREX::ERROR: invalid nucleus choice "<<endl;
    return 0;
  }
 
  Float_t calcrsec, mott, form_factor;

  /*Mott cross section for point-like scattering = 
    (alpha*Z* (hc/2pi)*cos(theta/2))^2 / 400E^2(sin(theta/2))^4 */
  
  Float_t pi = 3.1415926; 
  Float_t halfangle_rad = (angle/2)*(pi/180);
  
  Float_t sin4 = pow(sin(halfangle_rad),4);
  Float_t znuc;
  if (nuc == 0) {
    znuc = 82;
  } else {
    znuc = 6;
  }

  mott = pow(((znuc*0.197*cos(halfangle_rad))/137),2)/ (400*pow(energy,2)*sin4);

  // qsq comes from hamcKine 

  if (nuc == 0) {  // lead

    vector<Float_t>::const_iterator iterQsq = upper_bound(qsq_row.begin(), qsq_row.end(), qsq); //search for the first value of qsq which is larger or equal than actual   
    int indxQsq = iterQsq - qsq_row.begin(); 

    if (indxQsq <= 0 || indxQsq >= (Int_t)qsq_row.size()) return 0;
    
    Float_t qsq1, qsq2, form_factor1, form_factor2;
    qsq1 =  qsq_row.at(indxQsq);
    qsq2 = qsq_row.at(indxQsq-1);

    form_factor1 = ffsq_row.at(indxQsq);
    form_factor2 = ffsq_row.at(indxQsq-1);
    form_factor = Interpolate(qsq1, qsq2, qsq, form_factor1, form_factor2);

  } else { // carbon

    vector<Float_t>::const_iterator iterQsq = upper_bound(c12qsq_row.begin(), c12qsq_row.end(), qsq); //search for the first value of qsq which is larger or equal than actual   
    int indxQsq = iterQsq - c12qsq_row.begin(); 

    if (indxQsq <= 0 || indxQsq >= (Int_t)qsq_row.size()) return 0;

    Float_t qsq1, qsq2, form_factor1, form_factor2;
    qsq1 =  c12qsq_row.at(indxQsq);
    qsq2 = c12qsq_row.at(indxQsq-1);

    form_factor1 = c12ffsq_row.at(indxQsq);
    form_factor2 = c12ffsq_row.at(indxQsq-1);
    form_factor = Interpolate(qsq1, qsq2, qsq, form_factor1, form_factor2);

  }

  calcrsec = mott*form_factor; //result is in barn/seradians multiply by 1000 to compare with the value from Horowitchs table where it is milibars/stereadians 
  
  return calcrsec;
}

Float_t hamcPhyPREX::CalculateAsymmetry(Int_t nuc) {

// nuc = 0 --> lead
// nuc = 1 --> C12
// no other choice

  /*A(LR) = (G * Q^2 )/ (4pi*alpha*Sqrt(2))* [1-4sin(theta w)^2 - NeutronFF/ProtonFF ]
  GF = 1.6637 * 10^-5 GeV^-2 Fermi constant
  alpha = 1/137
  sin(theta w)^2 = 0.227 Weinberg Angle

  1 - 4*sin^2(theta_W) - N/Z = -1.44
  G / (4 pi alpha sqrt(2)) = 128.3 ppm / GeV^2  
  Asymmetry = - 128.3 * -1.44 * Q^2 = +184.74 Q^2 ppm
  But if N=Z, then Asymmetry = -128.3 * Q^2 ppm
  If lead, we'll use N=208, Z=82 and if C12, N=Z.
  
  */

  // qsq comes from hamcKine

  if (nuc != 0 && nuc != 1) {
    cout << "hamcPhyPREX::ERROR: incorrect nuc choice (0 or 1)"<<endl;
    return 0;
  }

  Float_t A0;
  if (nuc == 0) {  
    A0 = 184.74;  // lead
  } else {
    A0 = 128.3;   // C12
  }

  Float_t xasy = A0*qsq*pow(10.0,-6);
  asy0 = xasy;
  asy1 = xasy;  // no model dependence here.

  return xasy;
}

Float_t hamcPhyPREX::CalculateDrate(Float_t anum, Float_t tdens,Float_t tlen, Float_t crsec) {

  tlen = tlen*100;   // need cm
  //cout<<"tlen="<<tlen<<", tdens = "<<tdens<<", anum="<<anum<<", crsec = "<<crsec<<endl;
  Float_t avg_omega = 0.004671;
  Float_t xdrate = 6.25e12 * crsec * 0.602 * tlen * tdens * avg_omega / anum; //Hz/uA
  return xdrate;
}

Int_t hamcPhyPREX::LoadFormFactorTable(){

  /*Load Lead form factor
   first number is q, fifth is form factor squared 
   vector<Float_t> qsq_row, ffsq_row;  
   q = 0.197 * q;
   qsq = q*q;
  */
  
  FILE *fd;
  char strin[200]; 
  char* filename = "mefcal.pb208_1_25deg_fine.out"; 

  float q, ffsq, qsq;
  float ignore0, ignore1, ignore2, ignore3, ignore4, ignore5, ignore6, ignore7;
  
  fd = fopen(filename, "r");
  if (fd==NULL) {
    printf("ERROR: file %s does not exist \n", filename);
    printf("Bye Bye. \n");
    exit(0);
  }

  while(fgets(strin,1000,fd)!=NULL) {
    sscanf(strin, "%f %f %f %f %f %f %f %f %f %f", &ignore0, &q, &ignore1, &ignore2, &ffsq, &ignore3, &ignore4, &ignore5, &ignore6, &ignore7);

    q = 0.197*q; 
    qsq = pow(q, 2);
    qsq_row.push_back(qsq);
    ffsq_row.push_back(ffsq);
  }

  fclose(fd);

  return 1;
}


Int_t hamcPhyPREX::LoadHorowitzTable(vector<vector<Float_t> >& crsc_table, vector< vector<Float_t> >& asymmetry_table, Int_t stretch){

/*This function reads Horowitz's tables into two 2D vectors, one for crossection, another for asymmetry*/

  vector<float> crsc_row;  //temporary variables
  vector<float> asymmetry_row;

  Int_t debug = 0;

  energy_row.clear();
  asymmetry_row.clear();

  FILE *fd;
  float energy, angle, cross_section, asymmetry;
  float ignore1, ignore2, ignore3;
  char strin[100];

  /*Check wich filename*/
  char* filename;

  if (whichmodel == HORPB) {
     if (stretch==0) {
       filename="./PREX/horpb.dat";
     }  else {
       filename="./PREX/horpb1.dat";
     }
  }
  if (whichmodel == HORCA) {
     if (stretch==0) {
       filename="./PREX/ca48.dat";
     }  else {
       filename="./PREX/ca48s.dat";
     }
  }
  if (whichmodel == HORSN) {
     if (stretch==0) {
       filename="./PREX/sn120.dat";
     }  else {
       filename="./PREX/sn120s.dat";
     }
  }
  if (whichmodel == SI) filename = "./PREX/si.dat";
  if (whichmodel == NL3P06) filename = "./PREX/nl3p06.dat";
  if (whichmodel == SLY4) filename = "./PREX/sly4.dat";
  if (whichmodel == SIII) filename = "./PREX/siii.dat";
  if (whichmodel == FSU) filename = "./PREX/fsu.dat";
  if (whichmodel == NL3) filename = "./PREX/nl3.dat";
  if (whichmodel == NL3M05) filename = "./PREX/nl3m05.dat";

  cout << "PREX : using model "<<whichmodel<<"  stretch = ";
  cout << stretch << "    filename = "<< filename<<endl;

  fd=fopen(filename, "r");
  if (fd==NULL) {
    printf("hamcPhyPREX::ERROR: file %s does not exist \n", filename);
    printf("Bye bye !! \n");
    exit(0);
  }

  bool isFirst = true;
  bool didWriteAngle = false;
  while(fgets(strin,100,fd)!=NULL) {
    if (strstr(strin, "E=")!=NULL) { //Line for energy
      sscanf(strin, "E=%f", &energy);
      if (!isFirst) {
	crsc_table.push_back(crsc_row); /*if it is not the first energy
					  add row with data*/
	asymmetry_table.push_back(asymmetry_row);
	didWriteAngle = true;
	/* This is not first Energy value, therefore it has gone through angles and saved in the angle_row */
      } else {
	isFirst = false;	
      }
      energy_row.push_back(energy/1000); /*divide by 1000,because energies in Horowitch t tables are in MeV, but energy generated is in GeV  */	 
	crsc_row.clear(); //empty temporary variable
	asymmetry_row.clear();
    } else {
      // Format depends on the model
      if (whichmodel==HORPB || whichmodel==FSU) {
         sscanf(strin, "%f %f %f %f %f %f",
	     &angle, &cross_section, &ignore1, &asymmetry, &ignore2, &ignore3);
      } else {
        sscanf(strin, "%f %f %f ",
	       &angle, &cross_section, &asymmetry);
      }
      if (didWriteAngle==0) {
	angle_row.push_back(angle);	
      }
      crsc_row.push_back(cross_section);

      asymmetry_row.push_back(asymmetry);
      if (debug) {
       cout << "Loading  ... E = "<<energy<<"  "<<cross_section<<"  "<<crsc_row.size()<<"   "<<asymmetry<<"  "<<asymmetry_row.size()<<endl;
      }
    }
  }
  crsc_table.push_back(crsc_row);
  asymmetry_table.push_back(asymmetry_row);
  //  cout << "crsc_table size "<<crsc_table.size()<<"    asym tab size "<<asymmetry_table.size()<<endl;
  
  fclose(fd);

  return 1;
}

Int_t hamcPhyPREX::LoadC12FormFactorTable(){

  /*Load C12 form factor, needed for diamond foils.
   first number is q, fifth is form factor squared 
   vector<Float_t> qsq_row, ffsq_row;  
   q = 0.197 * q;
   qsq = q*q;
  */
  
  FILE *fd;
  char strin[200]; 
  char* filename = "mefcal.c12_1_25deg_fine.out"; 

  float q, ffsq, qsq;
  float ignore0, ignore1, ignore2, ignore3, ignore4, ignore5, ignore6, ignore7;
  
  fd = fopen(filename, "r");
  if (fd==NULL) {
    printf("ERROR: file %s does not exist \n", filename);
    printf("Bye Bye. \n");
    exit(0);
  }

  while(fgets(strin,1000,fd)!=NULL) {
    sscanf(strin, "%f %f %f %f %f %f %f %f %f %f", &ignore0, &q, &ignore1, &ignore2, &ffsq, &ignore3, &ignore4, &ignore5, &ignore6, &ignore7);

    q = 0.197*q; 
    qsq = pow(q, 2);
    c12qsq_row.push_back(qsq);
    c12ffsq_row.push_back(ffsq);
  }

  fclose(fd);

  return 1;
}
