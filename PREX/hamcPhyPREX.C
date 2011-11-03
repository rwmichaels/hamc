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
  didWriteAngle = false;
  num_models = 2;  // if =2 we're considering stretched vs unstretched R_N
  qsq = 0;
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
 
  //  PrintAsymTable();


  didinit = kTRUE;


  if (histo_test) {

    //    Float_t x1,x2;
    //    for (Int_t j=0; j<100; j++) {
    //       x1 = 0.8 + 0.4*((Float_t)j)/100;
    //       x2 = round(0.98*x1*20)/20;
    //       cout << "round test "<<j<<"   "<<x1<<"   "<<x2<<endl;
    //    }

    hpph1 = new TH1F("hpph1","Horowitz Lookup Crsec  vs  angle",2000,3,9);
    hpph2 = new TH1F("hpph2","Mott * FF (1.05 GeV)",2000,3,9);
    hpph3 = new TH1F("hpph3","Cross Section (0.2482 GeV)",2000,17,95);
    hpph4 = new TH1F("hpph4","Cross Section (0.502 GeV)",2000,12,42);
    hpph5 = new TH2F("hpph5","Percent diff vs angle",100,3,9,100,-10,2);
    hpph6 = new TH1F("hpph6","Asymmmetry  vs  angle",2000,3,9); 

    Float_t energy[4];
    energy[0] = 1.063;
    energy[1] = 1.063; 
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
        cout << "Energy "<<energy[iene]<<"     theta "<<theta_degr<<"     Crsec:   "<<crsec<<"          "<<crsec2<<"    "<<diff<<endl;

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

    Float_t theta_degr, theta_rad, ene, crsec2;

    for (Int_t itry=0; itry<20; itry++) {

      theta_degr = 5.0;
      Float_t E0nom = 1.05;
      ene = E0nom;

      if (itry==1) {
	theta_degr = 3.5;
        ene = E0nom;
      }
      if (itry==2) {
	theta_degr = 4.0;
        ene = E0nom;
      }
      if (itry==3) {
	theta_degr = 4.2;
        ene = E0nom;
      }
      if (itry==4) {
	theta_degr = 4.4;
        ene = E0nom;
      }
      if (itry==5) {
	theta_degr = 4.6;
        ene = E0nom;
      }
      if (itry==6) {
	theta_degr = 4.8;
        ene = E0nom;
      }
      if (itry==7) {
	theta_degr = 5.0;
        ene = E0nom;
      }
      if (itry==8) {
	theta_degr = 5.2;
        ene = E0nom;
      }
      if (itry==9) {
	theta_degr = 5.4;
        ene = E0nom;
      }
      if (itry==10) {
	theta_degr = 5.6;
        ene = E0nom;
      }
      if (itry==11) {
	theta_degr = 5.8;
        ene = E0nom;
      }
      if (itry==12) {
	theta_degr = 6.0;
        ene = E0nom;
      }
      if (itry==13) {
	theta_degr = 6.2;
        ene = E0nom;
      }
      if (itry==14) {
	theta_degr = 6.4;
        ene = E0nom;
      }
      if (itry==15) {
	theta_degr = 6.6;
        ene = E0nom;
      }
      if (itry==16) {
	theta_degr = 6.8;
        ene = E0nom;
      }
      if (itry==17) {
	theta_degr = 7.0;
        ene = E0nom;
      }
      if (itry==18) {
	theta_degr = 5.0;
        ene = 1.055;
      }
      if (itry==19) {
	theta_degr = 5.0;
        ene = 1.050;
      }

      theta_rad = 3.1415926*theta_degr/180.0;
      CrossSection(ene,theta_rad,0);
      Asymmetry(ene,theta_rad,0);
      Asymmetry(ene,theta_rad,1);

      Float_t frecoil = 1 + (ene/195.)*(1-TMath::Cos(theta_rad));
      Float_t eprime = ene/frecoil;
      qsq = 2*ene*eprime*(1-TMath::Cos(theta_rad));

      cout << "CHECK :  energy "<<ene<<"  theta "<<theta_degr<<endl;
      cout << "crsec  "<<crsec<<"    A0= "<<asy0<<"  A1= "<<asy1<<endl;
// Warning: CalculateCrossSection needs qsq, a variable member of this class.
      crsec2 = CalculateCrossSection(0, ene, theta_degr);  

      cout << "crsec2 = "<<crsec2<<endl;
      Float_t xdif = (crsec - crsec2)/crsec;
      cout << "Rel. diff  = "<<xdif<<"   ratio = "<<crsec/crsec2<<endl;

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
      Float_t daa, rate, asyavg, fom;
      Float_t daasum, ratesum, asysum, xnorm;
      Float_t domega = 0.0037;

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
    


  if (quick_fom) {  

    
    Float_t dcor[27];

    // May 2, 2010: these corrections are B.S.   Don't use !
    dcor[0] = 0.38;  dcor[1] = 0.33;  dcor[2]  = 0.20;  dcor[3] = 0.34;
    dcor[4] = 0.49;  dcor[5] = 0.78;  dcor[6]  = 1.03;  dcor[7] = 1.12;
    dcor[8] = 1.22;  dcor[9] = 1.31;  dcor[10]  = 1.35;  dcor[11] = 1.39;
    dcor[12] = 1.51;  dcor[13] = 1.61;  dcor[14]  = 1.73;  dcor[15] = 1.93;
    dcor[16] = 2.0;  dcor[17] = 2.16;  dcor[18]  = 2.20;  dcor[19] = 2.41;
    dcor[20] = 2.74;  dcor[21] = 4.2;  dcor[22]  = 5.6;  dcor[23] = 8.4;
    dcor[24] = 12.;  dcor[25] = 14.;  dcor[26]  = 16.;  


    Int_t use_MCcorr = 0;  // to use MC correcton or not
    Float_t xnfact = 0.484;

    Int_t itgt = 0;
    if (whichmodel == HORPB) itgt = 1;

    Float_t angcut = -2;     // degrees

  
    hd1 = new TH1F("hd1","Diff of crsec to lookup",200,-1,1);
    hfom1 = new TH1F("hfom1","sensitivity vs energy",200,0.65,1.25);
    hfom2 = new TH1F("hfom2","FOM (a.u.) vs energy",200,0.65,1.25);
    hfom3 = new TH1F("hfom3","FOM (a.u.) vs qsq",200,0.004,0.01015);
    hfom4 = new TH1F("hfom4","Raw dR/R(%) vs qsq",200,0.004,0.01015);
    hfom5 = new TH1F("hfom5","Total dR/R(%) vs qsq",200,0.004,0.01015);
    hfom6 = new TH1F("hfom6","1-arm Rate(MHz) vs qsq",200,0.004,0.01015);
    hfom7 = new TH1F("hfom7","Sensitivity(dA/A %) vs qsq",200,0.004,0.01015);
    hfom8 = new TH1F("hfom8","Phys Asym (p_e=1) vs qsq",200,0.004,0.01015);
    hfom9 = new TH1F("hfom9","Avg qsq vs E",200,0.9,1.2);
    hfom10 = new TH1F("hfom10","dA/A (%) vs qsq",200,0.004,0.01015);
    hfom11 = new TH1F("hfom11","Old Acceptance Function",600,3,8);
    hfom12 = new TH1F("hfom12","Correction factor to Acceptance",600,3,8);
    hfom13 = new TH1F("hfom13","Data-Corrected Acceptance Function",600,3,8);


    Float_t tdens = expt->target->GetMtlDensity(itgt);  // tgt density (g/cm^3)
    Float_t tlen = expt->target->GetMtlEffLen(itgt);  // tgt len (m)
    tlen = tlen*100;                        // need cm
    Float_t current = expt->event->beam->beam_current;  // microAmps (uA)
    current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec
    Float_t anum = expt->target->GetAscatt();    // atomic num.
    Float_t pol = 0.85;
    Float_t polerr = 0.015;
    Float_t domega = 0.0037;
    Float_t datarate = 4.4e8;

    cout << "Check FOM vs Energy "<<endl;
    cout << "Tgt  A = "<<anum<<"    tdens = "<<tdens<<"   tlen = "<<tlen<<endl;
    cout << "Current = "<<current<<endl;
    cout << "Polarization = "<<pol<<"    polerr "<<polerr<<endl;

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

      Float_t theta_degr, ene, xacc, xcrsec;
      Float_t daa, rate, asyavg, fom;
      Float_t daasum, ratesum, asysum, xnorm,xdiff;
      Float_t qsqsum, qsqavg, xcnt, asy_err, drr, sensi;
      Float_t blowup, drrtot;

      Float_t x1sum=0;  
      Float_t x2sum=0;

      for (Int_t iene=0; iene < 96; iene++) {
 
        ene = 0.72 + 0.005*((Float_t)iene); // GeV

        ratesum = 0;
        xnorm=0;
        asysum = 0;
        daasum = 0;
        qsqsum = 0;

        fclose(fd);
        fd=fopen("./PREX/accept.dat", "r");

        while(fgets(schar,100,fd)!=NULL) {
          sscanf(schar, "%f  %f", &theta_degr, &xacc);

	  //          if (theta_degr < angcut) continue;

	  // offset to simulate shift in Qsq
	  //  theta_degr += 0.25;

          theta_rad = 3.1415926*theta_degr/180.0;

          CrossSection(ene,theta_rad,0);  // uses lookup

          Float_t frecoil = 1 + (ene/195.)*(1-TMath::Cos(theta_rad));
          Float_t eprime = ene/frecoil;
          qsq = 2*ene*eprime*(1-TMath::Cos(theta_rad));
          xcrsec = CalculateCrossSection(0, ene, theta_degr);
          xdiff = xcrsec - crsec;
          hd1->Fill(xdiff); 

  // Lookup correction to MC
          Float_t xcor = 1;

          for (Int_t icc=0; icc<26; icc++) {
	    Float_t anglo = 0.068 + ((Float_t)(icc-1)) * 0.002;
            if (theta_rad > anglo && theta_rad < anglo+0.002) {
	      xcor = dcor[icc] + ((theta_rad-anglo)/0.002) * (dcor[icc+1]-dcor[icc]);
              xcor = xcor * xnfact;
              goto done1;
	    }
	  }
 done1:
          if (iene == 0) {
	    hfom11->Fill(theta_degr,xacc);
	    hfom12->Fill(theta_degr,xcor);
	    hfom13->Fill(theta_degr,xcor*xacc);
            x1sum += xacc;
            x2sum += xcor*xacc;
	  }

          daa = InterpAsym(ene, theta_rad);  // This sets asy0 and asy1 too

          Float_t xfact = xacc;
          if (use_MCcorr) xfact = xfact * xcor;

          rate = xfact * current * xcrsec * 0.602 * tlen * tdens * domega / anum;
         
          ratesum += rate;  
          daasum += daa * rate;
          asysum += asy0 * rate;
          qsqsum += qsq * rate;
          xnorm += xacc;

	  if (ene > 0.9 && ene <= 1.2) { 
	     if (use_MCcorr) {
	        cout << "dont use MCcorr for now";
                exit(0);
	     }
	     cout << "energy "<<ene<<"  theta "<<theta_degr<<endl;
             cout << "qsq "<<qsq<<endl;
	     cout << "accept  "<<xacc<<"  A= "<<asy0<<"  "<<asy1<<"  daa = "<<daa<<endl;
	     cout << "crsec  "<<crsec<<"   "<<xcrsec<<endl;
	     cout << "rate "<<rate<<"  "<<ratesum<<"  cnt = "<<xnorm<<endl;
             cout << "data rate "<<datarate<<endl;
	  }
	}
        if (ratesum > 0) {        
          asyavg = asysum / ratesum;  // does NOT have the p_e in it yet
          sensi = daasum / ratesum;
          rate = ratesum / xnorm;
          qsqavg = qsqsum / ratesum;
          fom = 1e9 * rate * (asyavg*asyavg) * (sensi*sensi);
// For 2 HRS and 30 days.  Use datarate (observed) instead of MC rate
          xcnt = 2.0 * datarate * 34. * 24. * 3600;
          asy_err = 0;   
          if (xcnt != 0) asy_err = 1.0 / TMath::Sqrt(xcnt);
          daa = asy_err / (pol * asyavg);  // stat. error.
          drr = 1;
          if (sensi != 0) drr = 0.01 * daa / sensi;
          blowup = 1;
          if (daa != 0) blowup = (TMath::Sqrt(polerr*polerr + daa*daa))/daa;
          drrtot = drr * blowup;
	} else {
          asyavg=0; sensi=0; rate=0; fom=0; drr=1;
	}
        cout << "E= "<<ene<<"  <qsq> "<<qsqavg<<"   <A> = "<<1e6*asyavg<<" ppm";
        cout << "  rate = "<<rate<<"  actuallly "<<datarate<<" Hz     sensi = "<<sensi<<endl;
        cout << "daa(stat) "<<daa<<"    dR/R   "<<drr<<"   dR/R tot =  "<<drrtot<<endl;
        hfom1->Fill(ene,sensi);
        hfom2->Fill(ene,fom);
        hfom3->Fill(qsqavg,fom);
        hfom4->Fill(qsqavg,drr);
        hfom5->Fill(qsqavg,drrtot);
        hfom6->Fill(qsqavg,1e-6*rate);
        hfom7->Fill(qsqavg,-100*sensi);
        hfom8->Fill(qsqavg,asyavg);
        hfom9->Fill(ene,qsqavg);
        hfom10->Fill(qsqavg,daa);
      }

      cout << "sums xxx "<<x1sum<<"  "<<x2sum<<endl;
    }


  }  // quick_fom
   

  if (power_integ) {  

    // To integrate the power from anglo_deg to anghi_deg


    Float_t theta_degr, dtheta, steptheta, ene, xcrsec;
    Float_t rate, domega, xpts, ratesum;
    Float_t radlen, prob, xfact, xinv, eloss; 
    
    Int_t Npts = 1000;
    xpts = ((Float_t)Npts);

    Float_t tdens = 11.35;          // g/cm^3
    Float_t tlen = 0.05;            // cm
    Float_t current = 100;          // uA
    current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec
    Float_t anum = 208;             // atomic num.
    Float_t pi = 3.1415926;

    Float_t anglo_deg = 10;  // degrees
    Float_t anghi_deg = 30;

    Int_t use_brehms = 1;  // DONT TOUCH

    cout << "Check Power in solid angle from "<<endl;
    cout <<  anglo_deg << "   to  "<<anghi_deg<<"   degrees "<<endl;
    cout << "At energy "<<ene<<endl;
    cout << "Tgt  A = "<<anum<<"    tdens = "<<tdens<<"   tlen = "<<tlen<<endl;
    cout << "Current = "<<current<<endl;

    steptheta = (anghi_deg-anglo_deg)/xpts;
    dtheta = (pi/180.)*steptheta;

    ratesum = 0;

    for (Int_t iang = 0; iang<Npts; iang++) {

      ene = 1.063;  // GeV

      theta_degr = anglo_deg + steptheta * ((Float_t)iang);

      theta_rad = (pi/180.)*theta_degr;
  
      if (use_brehms) {
           prob = 0.0000000001+gRandom->Rndm();
           if (prob>1.0) prob = 1;
           radlen = 0.05; // 1/2 target
           xinv = 1./radlen;
           xfact = TMath::Exp(xinv*TMath::Log(prob));
           eloss = ene*xfact;
	   // test
           if (eloss > 0.8) eloss = 0.8;
           ene = ene - eloss;
           cout << "Brehms "<<xfact<<"  "<<eloss<<"  "<<ene<<endl;
      }

      Float_t frecoil = 1 + (ene/195.)*(1-TMath::Cos(theta_rad));
      Float_t eprime = ene/frecoil;
      qsq = 2*ene*eprime*(1-TMath::Cos(theta_rad));
      xcrsec = CalculateCrossSection(0, ene, theta_degr);
      domega = 2 * pi * TMath::Sin(theta_rad) * dtheta;

      rate = current * xcrsec * 0.602 * tlen * tdens * domega / anum;

      // Weight the rate by energy 
      rate = rate * ene/1.063;

      ratesum += rate;  

      cout << "energy "<<ene<<"  theta "<<theta_degr<<endl;
      cout << "crsec  "<<xcrsec<<"   domega "<<domega<<endl;
      cout << "rate "<<rate<<"  "<<ratesum<<endl;

    }

    cout << endl << "Total rate "<<ratesum<<"  for current "<<current/6.25e12<<"   uA"<<endl;
    cout << "Power = "<<ratesum*ene*1.6e-10<<"  Watts "<<endl;

    exit(0);


  }  // power_integ
   

  if (check_ms_1D) {  

    // Crude check of acceptance in 1D

    Int_t ldebug = 0;

    cout << "Checking MS in 1D"<<endl;

    Float_t ene, xcrsec;
    Float_t frecoil, eprime, qsq_ms;
    Float_t relrate, domega;
    Float_t theta_deg, theta_rad, dtheta_rad, theta_ms_rad;

    Int_t inaccept;

    Int_t Nang=100000;
    Float_t anglo_deg, anghi_deg;   // degrees
    Float_t anglo = 1;            // lowest angle
    Float_t anghi = 10.0;         // highest 
    Float_t angbite, xang;
    Float_t angle_nom_deg = 5.0;   // HRS nominal angle
    Float_t angle_nom_rad;
    Float_t sth, cth, tth;
    Float_t trk_sth, trk_cth, trk_tth;
    Float_t xtrk, ytrk;
    Float_t xctr,yctr,x1,y1,x2,y2;

    Float_t zlen = 2.5;   // meters (dist of target to collimator)
    Float_t zcol = 0.06;   // size of collimator (meters)

    Float_t pi = 3.1415926;

    histang0 = new TH1F("histang0","Accepted Angles in 1D (degr) ",100,0.2,10);
    histang1 = new TH1F("histang1","Mult. Scatt. (rad) ",500,-0.05,0.05);
    histang2 = new TH1F("histang2","Accept Fcn (rel. to central ang.) no MS",100,-0.1,0.1);
    histang3 = new TH1F("histang3","Accept Fcn (rel. to central ang.) with MS",100,-0.1,0.1);
    histang4 = new TH1F("histang4","Accept Fcn (rel. to central ang.) no MS, weighted",100,-0.1,0.1);
    histang5 = new TH1F("histang5","Accept Fcn (rel. to central ang.) with MS, weighted",100,-0.1,0.1);
    histang6 = new TH1F("histang6","qsq in accept",100,0,0.02);
    histang7 = new TH1F("histang7","qsq_obs in accept",100,0,0.02);

    angle_nom_rad = (pi/180.)* angle_nom_deg;
    sth = TMath::Sin(angle_nom_rad);
    cth = TMath::Cos(angle_nom_rad);
    tth = sth/cth;

    xctr = zlen * sth;
    yctr = zlen * cth;
    x1 = xctr - zcol/cth;
    x2 = xctr + zcol/cth;
    y1 = yctr + -1*tth*(x1-xctr);  // note, slope is negative.
    y2 = yctr + -1*tth*(x2-xctr);  //  "  "

    xang = (Float_t) Nang;
    angbite = (anghi - anglo) / xang;

    ene = 1.063;  // GeV

    for (Int_t iang = 0; iang < Nang; iang++) {

      anglo_deg = anglo + angbite * ((Float_t)iang);

      anghi_deg = anglo_deg + angbite;

      theta_deg = (0.5 * (anglo_deg + anghi_deg));
      theta_rad = (pi/180.)* theta_deg;
      dtheta_rad = (pi/180.)* angbite;

      frecoil = 1 + (ene/195.)*(1-TMath::Cos(theta_rad));
      eprime = ene/frecoil;
      qsq = 2*ene*eprime*(1-TMath::Cos(theta_rad));  // warning: qsq must be defined
      xcrsec = CalculateCrossSection(0, ene, theta_deg);
      domega = 2 * pi * TMath::Sin(theta_rad) * dtheta_rad;
      relrate = 10000 * xcrsec * domega;

// Multiple scattering in lead

      Float_t radlen = 0.1;
      Float_t theta_sigma = (0.0136/ene) * TMath::Sqrt(radlen);
      Float_t deltaMS, prob, prob2, xsign;

      prob = gRandom->Rndm();

      if (prob < 0.02) {  // flat tail
  	  xsign = 1;
          prob2 = gRandom->Rndm();
          if (prob2 < 0.5) xsign = -1;
          deltaMS = xsign * 10 * theta_sigma * gRandom->Rndm(); 
      } else {
          deltaMS = theta_sigma * gRandom->Gaus();
      }

      theta_ms_rad = theta_rad + deltaMS;
      qsq_ms = 2*ene*eprime*(1-TMath::Cos(theta_ms_rad));  

      histang1->Fill(deltaMS);

      if (theta_ms_rad <= 0) continue;

      trk_sth = TMath::Sin(theta_ms_rad);
      trk_cth = TMath::Cos(theta_ms_rad);
      trk_tth = trk_sth/trk_cth;

      xtrk = (1./((1./trk_tth) + tth)) * zlen * (cth + tth*sth);
      ytrk = xtrk/trk_tth;

      inaccept = 0;

      if ( xtrk > x1 && xtrk < x2 && ytrk < y1 && ytrk > y2) inaccept = 1;

      if (ldebug && (theta_deg > 4 && theta_deg < 6)) {
        cout << endl<<"angle of track "<<theta_deg<<endl;
        cout << "collim "<<zlen<<"  "<<zcol<<endl;
        cout << "Angle factors (nom) "<<sth<<"  "<<cth<<"  "<<tth<<endl;
        cout << "Angle factors (trk) "<<trk_sth<<"  "<<trk_cth<<"  "<<trk_tth<<endl; 
	cout << "x1,y1  "<<x1<<"   "<<y1<<endl;
	cout << "x2,y2  "<<x2<<"   "<<y2<<endl;
        cout << "xtrk, ytrk "<<xtrk<<"  "<<ytrk<<endl;
        cout << "relrate = "<<relrate<<"  "<<xcrsec<<"  "<<domega<<endl;
        if (inaccept) cout << "--------------------- IN ACCEPTANCE "<<endl;
      }

      if (inaccept) {

          histang0->Fill(theta_deg);

   	  histang2->Fill(theta_rad - angle_nom_rad);  // no MS
       
	  histang3->Fill(theta_ms_rad - angle_nom_rad);  // with MS

	  histang4->Fill(theta_rad - angle_nom_rad, relrate);  // no MS, weight by rate
       
	  histang5->Fill(theta_ms_rad - angle_nom_rad, relrate);  // with MS, weight by rate

          histang6->Fill(qsq, relrate);

          histang7->Fill(qsq_ms, relrate);

      }
    }
 
  }  // check_ms_1D


  if (neutron_power) {  

    // To integrate the power differentially with angle

    Int_t ldebug = 0;

    Int_t use_multsc = 1;
    Int_t use_brehms = 1;
    
    Float_t theta_degr, dtheta_rad, steptheta_degr;
    Float_t ene0, ene, xcrsec, xnorm;
    Float_t frecoil, eprime;
    Float_t rate, domega, xang, xsim, xcnt, xcnt2, xpts;
    Float_t ratesim, ratesimdump, ratewhat;
    Float_t ratep0, ratep1, ratep2, ratep3, ratep4;
    Float_t radlen, prob, ptot, xfact, xinv, eloss; 
    Float_t xt, xpt, yt, ypt, xt0, xpt0, yt0, ypt0;
    Float_t phi, rrad;
    Float_t zdrift1, zdrift2, zdrift3;
    Float_t radcut1, radcut2, radcut3, radcut4, radcut5;
    Float_t xptcut1, yptcut1, xptcut2, yptcut2; 
    Float_t xptcut3, yptcut3, xptcut4, yptcut4; 
    Int_t inacc;

    histphi = new TH1F("histphi","Phi chk",100,-0.2,6.5);
    histrast = new TH2F("histrast","Rast chk",100,-1,1,100,-1,1);
    histene  = new TH1F("histene","Energy ",100,0.05,1.2);
    histprob = new TH1F("histprob","prob vs angle",100,0.0,1.0);
    histinacc = new TH1F("histinacc","In acceptance",10,-0.5,2.0);

    histz1a = new TH1F("histz1a","Radius (cm) at entrance to SEPTUM",100,-0.2,100.0);
    histz1b = new TH1F("histz1b","X angle at Z1",100,-0.17,0.17);
    histz1c = new TH1F("histz1c","Y angle at Z1",100,-0.17,0.17);
    histz1d = new TH1F("histz1d","(inacc) Radius (cm) at Z1",100,-0.2,100.0);
    histz1e = new TH1F("histz1e","(inacc) X angle at Z1",100,-0.17,0.17);
    histz1f = new TH1F("histz1f","(inacc) Y angle at Z1",100,-0.17,0.17);

    histz2a = new TH1F("histz2a","Radius (cm) at ",100,-0.2,100.0);
    histz2b = new TH1F("histz2b","X angle at Z2",100,-0.17,0.17);
    histz2c = new TH1F("histz2c","Y angle at Z2",100,-0.17,0.17);
    histz2d = new TH1F("histz2d","(inacc) Radius (cm) at Z2",100,-0.2,100.0);
    histz2e = new TH1F("histz2e","(inacc) X angle at Z2",100,-0.17,0.17);
    histz2f = new TH1F("histz2f","(inacc) Y angle at Z2",100,-0.17,0.17);

    histz3a = new TH1F("histz3a","Radius (cm) after SEPTUM and DRIFT",100,-0.2,100.0);
    histz3b = new TH1F("histz3b","(inacc) Radius (cm) at Z3",100,-0.2,100.0);

    histz4a = new TH1F("histz4a","Radius (cm) after COMPENSATING QUAD",100,-0.2,100.0);
    histz4b = new TH1F("histz4b","(inacc) Radius (cm) at Z4",100,-0.2,100.0);

    histz5a = new TH1F("histz5a","Radius (cm) at BEAM DUMP ",100,-0.2,150.0);
    histz5b = new TH1F("histz5b","Radius (cm) ad dump (cut) ",100,-0.2,150.0);

    Int_t Nang = 100;   // number of angle bins
    Int_t Npts = 200;   // number of points within a bin to integ

    Int_t Nsim = 100;  // number of points per angle bin to simulate (Brehms, and later phi and raster)

    quad1 = new hamcQuad(800, 1, 87);
    quad1->SetFocusY();
    quad2 = new hamcQuad(200, 1, 87);
    quad2->SetFocusX();

    cout << "Quad 1 setup "<<endl;
    quad1->Print();
    cout << "\n\nQuad 2 setup "<<endl;
    quad2->Print();

    // all units cm
    // to septum     between quads    to dump
    zdrift1 = 150;   zdrift2 = 60;   zdrift3 = 2100;

// Transport model is D Q D Q D to dump.  
// radcut after each each element (units: cm)                       
    radcut1 = 12;  radcut2 = 20;   radcut3 = 30;  radcut4 = 40;  radcut5 = 80;

    xptcut1 = 0.08; 
    xptcut2 = 0.08; 
    xptcut3 = 0.08; 
    xptcut4 = 0.08; 
    yptcut1 = 0.08; 
    yptcut2 = 0.08; 
    yptcut3 = 0.08; 
    yptcut4 = 0.08; 

    xang = ((Float_t)Nang);
    xsim = ((Float_t)Nsim);
    xpts = ((Float_t)Npts);

    Float_t tdens = 11.35;          // g/cm^3
    Float_t tlen = 0.05;            // cm
    Float_t current = 100;          // uA
    current = current * 6.25e12;    // 100 uA = 6.25e14 e- / sec
    Float_t anum = 208;             // atomic num.
    Float_t pi = 3.1415926;

    Float_t anglo_deg, anghi_deg;   // degrees
    Float_t anglolo = 0.002;        // lowest we'll go
    Float_t angbite;                // angle bite
    Float_t anghihi = 10.0;         // highest we'd go
    Float_t angp0_lo;               // lowest where we believe rate
    Float_t angp1_lo, angp1_hi;     // limits for PREX I (2010)
    Float_t angp2_lo, angp2_hi;     // poss. limits for PREX II 


    angp0_lo = anglolo;  // essentially no cut
 

    angp1_lo = 1.1;
    angp1_hi = anghihi;

    angp2_lo = 3.0;  // ???
    angp2_hi = anghihi;

    angbite = (anghihi - anglolo) / xang;

    ene0 = 1.063;  // GeV

// Neutron power

    cout << "Power integral for elastic scattering"<<endl;
    cout << "At energy "<<ene0<<endl;
    cout << "Tgt  A = "<<anum<<"    tdens = "<<tdens<<"   tlen = "<<tlen<<endl;
    cout << "Current = "<<current<<"  e-/sec =   "<<current/6.25e12<<"   uA"<<endl;

    ratep0 = 0;
    ratep1 = 0;
    ratep2 = 0;
    ratep3 = 0;
    ratep4 = 0;
    ptot = 0;

    for (Int_t iang = 0; iang < Nang; iang++) {

      anglo_deg = anglolo + angbite * ((Float_t)iang);

      anghi_deg = anglo_deg + angbite;

      if (anghi_deg > anghihi) goto done101;  // done

      if ((iang%10)==0) cout <<  anglo_deg << "   to  "<<anghi_deg<<"   degrees "<<endl;

      steptheta_degr = (anghi_deg-anglo_deg)/xpts;
      dtheta_rad = (pi/180.)*steptheta_degr;


// Integrate of over fine bins in angle within the bite "angbite".

      for (Int_t ibin = 0; ibin<Npts; ibin++) {

        theta_degr = anglo_deg + steptheta_degr * ((Float_t)ibin);

        theta_rad = (pi/180.)*theta_degr;

        ratesim = 0;
        ratesimdump = 0;
        ratewhat = 0;

// Within each scattering angle increment, 
// simulate by wiggling:  azimuth, raster, and Brehmstrahlung Eloss.

        for (Int_t ievt = 0; ievt<Nsim; ievt++) {

           ene = ene0;

  	   inacc = 1;

// Need a phi (azimuth)

           phi = 2*pi*gRandom->Rndm();

           histphi->Fill(phi);

           xpt = TMath::Tan(theta_rad) * TMath::Cos(phi);
           ypt = TMath::Tan(theta_rad) * TMath::Sin(phi);

// Need a raster (+/- 0.2 cm)

          xt  = -0.2 + 0.4*gRandom->Rndm();
          yt  = -0.2 + 0.4*gRandom->Rndm();

          histrast->Fill(xt,yt);
 
// Brehmstrahlung energy loss


          if (use_brehms) {

            prob = 0.0000000001+gRandom->Rndm();
            if (prob>1.0) prob = 1;
            radlen = 0.05; // 1/2 target
            xinv = 1./radlen;
            xfact = TMath::Exp(xinv*TMath::Log(prob));
            eloss = ene0*xfact;
        // truncate the Eloss
            if (eloss > 0.99*ene0) eloss = 0.99*ene0;
            ene = ene0 - eloss;
            if (ldebug) cout << "Brehms "<<xfact<<"  "<<eloss<<"  "<<ene<<endl;

	  }

          histene->Fill(ene);

          frecoil = 1 + (ene/195.)*(1-TMath::Cos(theta_rad));
          eprime = ene/frecoil;
          qsq = 2*ene*eprime*(1-TMath::Cos(theta_rad));
          xcrsec = CalculateCrossSection(0, ene, theta_degr);
          domega = 2 * pi * TMath::Sin(theta_rad) * dtheta_rad;
          prob = xcrsec * 0.602 * tlen * tdens * domega / anum;  
          ptot = ptot + prob;

          histprob->Fill(theta_degr, prob);

          rate = current * prob;
 
// Weight the rate by energy (power ~ energy)

          rate = rate * ene/1.063;

          ratesim += rate;  

          if (ldebug) {
            cout << "energy "<<ene<<"  theta "<<theta_degr<<endl;
            cout << "crsec  "<<theta_degr<<"  "<<xcrsec<<"   domega "<<domega<<"  "<<prob<<endl;
            cout << "rate in bin "<<rate<<endl;
            cout << "Initial trans " << xt << "  "<<yt<<"   "<<xpt<<"  "<<ypt<<endl;
	  }

// Multiple scattering in lead

          if (use_multsc) {

             Float_t radlen = 0.1;
             Float_t theta_sigma = (0.0136/ene) * TMath::Sqrt(radlen);

             Float_t dtheta1, dtheta2, prob;

             prob = gRandom->Rndm();
             if (prob < 0.02) {
                 dtheta1 = 5 * theta_sigma * gRandom->Rndm();  // flat tail
            } else {
                 dtheta1 = theta_sigma * gRandom->Gaus();
            }

            xpt = xpt + dtheta1;
      
            prob = gRandom->Rndm();
            if (prob < 0.02) {
               dtheta2 = 5 * theta_sigma * gRandom->Rndm();  // flat tail
            } else {
               dtheta2 = theta_sigma * gRandom->Gaus();
            }

            ypt = ypt + dtheta2;

	  }

 // Transport in x and y. D Q D Q D to dump.  Check aperture in each place. 

          xt  =  xt + xpt * zdrift1;
          yt  =  yt + ypt * zdrift1;

          rrad = TMath::Sqrt(xt*xt + yt*yt);

          if (ldebug) {
            cout << "z1 "<<zdrift1<<endl;
            cout << "Trans at z1 " << xt << "  "<<yt<<"   "<<xpt<<"  "<<ypt<<endl;
            cout << "radius at z1  "<<rrad << endl;
	  }

          if (rrad > radcut1 || xpt > xptcut1 || -1*xpt > xptcut1 ||
	       ypt > yptcut1 || -1*ypt > yptcut1 ) {
   	           inacc = 0;               
   	  }

          histz1a->Fill(rrad,prob);
          histz1b->Fill(xpt);
          histz1c->Fill(ypt);
          if (inacc==0) {
	    histz1d->Fill(rrad,prob);
	    histz1e->Fill(xpt);
	    histz1f->Fill(ypt);
	  }

          xt0 = xt;  xpt0 = xpt;
          xt  =  quad1->GetMatrix(ene,0) * xt0 + quad1->GetMatrix(ene,1) * xpt0;
          xpt =  quad1->GetMatrix(ene,2) * xt0 + quad1->GetMatrix(ene,3) * xpt0;

          yt0 = yt;  ypt0 = ypt;
          yt  =  quad1->GetMatrix(ene,4) * yt0 + quad1->GetMatrix(ene,5) * ypt0;
          ypt =  quad1->GetMatrix(ene,6) * yt0 + quad1->GetMatrix(ene,7) * ypt0;

          rrad = TMath::Sqrt(xt*xt + yt*yt);

          if (ldebug) {
            cout << "Quad 1  m.e. "<<quad1->GetMatrix(ene,0)<<"  "<<quad1->GetMatrix(ene,1);
            cout << "   "<<quad1->GetMatrix(ene,2)<<"  "<<quad1->GetMatrix(ene,3)<<endl;
            cout << "   "<<quad1->GetMatrix(ene,4)<<"  "<<quad1->GetMatrix(ene,5);
            cout << "   "<<quad1->GetMatrix(ene,6)<<"  "<<quad1->GetMatrix(ene,7)<<endl;
            cout << "Trans at z1 " << xt << "  "<<yt<<"   "<<xpt<<"  "<<ypt<<endl;
            cout << "radius at z1  "<<rrad << endl;
	  }

 
          if (rrad > radcut2) {
   	      inacc = 0;               
   	  }

          histz2a->Fill(rrad,prob);
          histz2b->Fill(xpt);
          histz2c->Fill(ypt);
          if (inacc==0) {
  	     histz2d->Fill(rrad,prob);
             histz2e->Fill(xpt);
             histz2f->Fill(ypt);
	  }

          xt  =  xt + xpt * zdrift2;
          yt  =  yt + ypt * zdrift2;

          rrad = TMath::Sqrt(xt*xt + yt*yt);

          if (rrad > radcut3 ) inacc = 0;               
   
          histz3a->Fill(rrad,prob);
          if (inacc==0) histz3b->Fill(rrad,prob);

          xt0 = xt;  xpt0 = xpt;
          xt  =  quad2->GetMatrix(ene,0) * xt0 + quad2->GetMatrix(ene,1) * xpt0;
          xpt =  quad2->GetMatrix(ene,2) * xt0 + quad2->GetMatrix(ene,3) * xpt0;

          yt0 = yt;  ypt0 = ypt;
          yt  =  quad2->GetMatrix(ene,4) * yt0 + quad2->GetMatrix(ene,5) * ypt0;
          ypt =  quad2->GetMatrix(ene,6) * yt0 + quad2->GetMatrix(ene,7) * ypt0;

          rrad = TMath::Sqrt(xt*xt + yt*yt);
 
          if (rrad > radcut4) inacc = 0;               

          histz4a->Fill(rrad,prob);
          if (rrad > radcut4) {
	    histz4b->Fill(rrad,prob);
          }

          xt  =  xt + xpt * zdrift3;  // transport to beam dump
          yt  =  yt + ypt * zdrift3;

          rrad = TMath::Sqrt(xt*xt + yt*yt);
          if (rrad > radcut5 ) inacc = 0;               

          histz5a->Fill(rrad,prob);
          if (inacc == 0 ) {
             histz5b->Fill(rrad,prob);
             ratesimdump = ratesimdump + rate;
	  } else {
 	     ratewhat = ratewhat + rate;
	  }

          histinacc->Fill(inacc);

	}

        rate = ratesim/xsim;

        ratep0 = ratep0 + rate;
 
        if (theta_degr > angp1_lo && theta_degr < angp1_hi) ratep1 = ratep1 + rate;

        if (theta_degr > angp2_lo && theta_degr < angp2_hi) ratep2 = ratep2 + rate;


        ratep3 = ratep3 + ratesimdump/xsim; 

        ratep4 = ratep4 + ratewhat/xsim; 

      }  // loop over i (angle bins)

      rate = ratesim/xsim;

    }

 done101:

    cout << endl << "Total rate "<<ratep0<<endl;
    cout << "Total prob "<<ptot<<endl;
    cout << "Power = "<<ratep0*ene*1.6e-10<<"  Watts "<<endl;
    cout << "Rate PREX I  =   "<<ratep1<<"  Power = "<<ratep1*ene*1.6e-10<<"  Watts "<<endl;
    cout << "Rate PREX II  =   "<<ratep2<<"  Power = "<<ratep2*ene*1.6e-10<<"  Watts "<<endl;
    cout << "Rate PREX II trans  =   "<<ratep3<<"  Power = "<<ratep3*ene*1.6e-10<<"  Watts "<<endl;
    cout << "Leftover "<<ratep4<<endl;

 
  }  // neutron_power
   

  return 1;
}


Float_t hamcPhyPREX::InterpAsym(Float_t ene, Float_t theta_rad) {

   Int_t ldebug = 0;

   Float_t daa = 1;

   Float_t asy0_lo, asy0_hi, asy1_lo, asy1_hi;

   Float_t eremain = 0.05 * ((Int_t)(ene/0.05));

   Asymmetry(eremain,theta_rad,0);
   Asymmetry(eremain,theta_rad,1);

   asy0_lo = asy0;
   asy1_lo = asy1;

   Asymmetry(eremain+0.05,theta_rad,0);
   Asymmetry(eremain+0.05,theta_rad,1);

   asy0_hi = asy0;
   asy1_hi = asy1;

   asy0 = asy0_lo + (asy0_hi-asy0_lo)*((ene-eremain)/0.05);
   asy1 = asy1_lo + (asy1_hi-asy1_lo)*((ene-eremain)/0.05);

   if (asy0 != 0) daa = (asy1 - asy0)/asy0;

   if (ldebug) {
     cout << "InterpAsym  "<<ene<<"  "<<theta_rad<<"   "<<eremain<<"  "<<endl;
     cout << asy0_lo << "   "<< asy0_hi << "   "<< asy1_lo << "   "<< asy1_hi;
     cout << "   "<< asy0 << "   "<<asy1<<endl;
   }

   return daa;

}


Int_t hamcPhyPREX::Generate(hamcExpt *expt) {

   Int_t ldebug=0;

   Float_t energy = kine->energy;    // GeV
   Float_t theta = kine->theta;      // radians
   qsq = kine->qsq;                  // GeV^2

   Float_t anum = expt->target->GetAscatt();

   Int_t mtl_idx = expt->target->GetMtlIndex();
   Float_t tdens = expt->target->GetMtlDensity(mtl_idx);  // tgt density (g/cm^3)
   Float_t tlen = expt->target->GetMtlLen(mtl_idx);  // tgt len (m)
   Float_t tefflen = expt->target->GetMtlEffLen(mtl_idx);  // eff. tgt len (m)

   if (ldebug) cout << "anum = "<<anum<<endl;

   if (anum == 12) {
     crsec = CalculateCrossSection(1, energy, theta*180/PI);
     asymmetry = CalculateAsymmetry(1);
     drate = CalculateDrate(anum, tdens, tlen, crsec);
     if (ldebug) {
       cout << endl<< "C12:  energy "<<energy<<"  theta "<<theta<<"  qsq "<<qsq<<endl;
       cout<< "C12 crsec = "<<crsec<<"   A = "<<asymmetry<<endl<<endl;
     }
     return OK;
     //     cout << "C12 "<<1e6*crsec<<endl;
   }

// Compute the cross section for PREX
   CrossSection(energy, theta, 0);
//   cout << "Pb  "<<1e6*crsec<<endl<<endl;

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

  if (debug) {
     cout << "\n\nenergy "<<energy<<"        angle "<<angle<<endl;
     cout << "Indices  "<<indxEnergy<<"   "<<indxAngle<<endl;
     cout << energy_row[indxEnergy] << "  " << crsc_tables[0][indxEnergy][0] << "    " << crsc_tables[0][indxEnergy][4] << endl;
  }

  Float_t crsc1, crsc2;
  Float_t crsca1, crsca2;
  Float_t e1, e2;

  //find the angles above and below the actual angle; used for interpolation
  Float_t angle_upper =angle_row[indxAngle+1];
  Float_t angle_lower = angle_row[indxAngle];

  if (debug) {
    cout << "angle lims "<<indxAngle<<"  "<<angle_lower<<"  "<<angle_upper<<endl;
  }

  crsc1 = crsc_tables[stretch][indxEnergy][indxAngle]; /*get the cross section value for angle value larger than the actual*/
  crsc2 = crsc_tables[stretch][indxEnergy][indxAngle+1]; //get the cross section value for nearby angle */
  crsca1 = Interpolate(angle_lower, angle_upper, angle, crsc1, crsc2)/1000;  

  if (debug) {
    cout << "crsec1, 2  "<<indxEnergy<<"  "<<crsc1<<"  "<<crsc2<<endl;
  }

  crsca2 = -1;

  if (indxEnergy < energy_row.size()-2) {
     crsc1 = crsc_tables[stretch][indxEnergy+1][indxAngle]; 
     crsc2 = crsc_tables[stretch][indxEnergy+1][indxAngle+1]; 
     crsca2 = Interpolate(angle_lower, angle_upper, angle, crsc1, crsc2)/1000;  
  }

  crsec = crsca1;

  if (crsca2 > -1) {

    e1 = energy_row[indxEnergy];
    e2 = energy_row[indxEnergy+1];     
    crsec = Interpolate(e1, e2, energy, crsca1, crsca2);  
  }

  if (debug) {
     cout << "\n LEAD Cross section : "<<endl;
     cout << "stretch "<<stretch<<"   indices "<<indxEnergy<<"   "<<indxAngle+1<<endl;
     cout << "energy " << energy << " GeV " << "angle " << angle << endl;
     cout << "e1 "<<e1<<"    e2 "<<e2<<endl;
     cout << "crsecs =  "<<crsca1<<"  "<<crsca2<<"  "<<crsec<<endl;
     
  }

  Float_t calccrsec = CalculateCrossSection(0, energy, angle);
  if (debug) {
       cout <<"calculated lead cross section " << calccrsec <<endl;
       cout <<" fractional diff "<< (calccrsec - crsec)/crsec<<endl;
  }

  return OK;
}


Int_t hamcPhyPREX::Asymmetry(Float_t energy, Float_t angle_rad, Int_t stretch) {

  if (!didinit) {
    cout << "ERROR: hamcPhyPREX:: Not initialized."<<endl;
       return -1;
  }
  // lookup the asymmetry for this E,theta

  Int_t ldebug = 0;

  Float_t angle = angle_rad * 180 / PI;

  asymmetry = 0;

  Int_t indxAngle = FindAngleIndex(angle);
  Int_t indxEnergy = FindEnergyIndex(energy);

  if (indxAngle <= 0 || indxEnergy <= 0) return -1;

  Float_t asymmetry1, asymmetry2, asyE1, asyE2, e1, e2;
  Float_t angle_upper = angle_row[indxAngle+1] ;
  Float_t angle_lower = angle_row[indxAngle];

  asymmetry1 = asymmetry_tables[stretch][indxEnergy][indxAngle];
  asymmetry2 = asymmetry_tables[stretch][indxEnergy][indxAngle+1];
  asyE1 = Interpolate(angle_lower, angle_upper, angle, asymmetry1, asymmetry2);

  if (ldebug) {
    cout << "\n\nEnergy  "<<energy<<"   angle "<<angle<<endl;
    cout << "indices "<<stretch<<"  "<<indxEnergy<<"  "<<indxAngle<<"  "<<indxAngle+1<<endl;
    cout << "asyE1 inputs "<<asymmetry1<<"  "<<asymmetry2<<"  "<<angle_lower<<"  "<<angle_upper<<"  "<<asyE1<<endl;
  }

  asyE2 = -1;

  if (indxEnergy < energy_row.size()-2) {
     asymmetry1 = asymmetry_tables[stretch][indxEnergy+1][indxAngle];
     asymmetry2 = asymmetry_tables[stretch][indxEnergy+1][indxAngle+1];
     asyE2 = Interpolate(angle_lower, angle_upper, angle, asymmetry1, asymmetry2);
  }

  if(ldebug) 
     cout << "asyE2 inputs "<<asymmetry1<<"  "<<asymmetry2<<"  "<<angle_lower<<"  "<<angle_upper<<"  "<<asyE2<<endl;

  asymmetry = asyE1;

  if (asyE2 > -1) {
    e1 = energy_row[indxEnergy];
    e2 = energy_row[indxEnergy+1];     
    asymmetry = Interpolate(e1, e2, energy, asyE1, asyE2);  
  }

  if (stretch == 0) asy0 = asymmetry;
  if (stretch == 1) asy1 = asymmetry;

  int debug=0;

  if (debug) {
    cout<<"angle  "<<angle<<"   "<<angle_lower<<"   "<<angle_upper<<endl;
    cout << "Indices "<<indxAngle<<"  "<<indxEnergy<<"   "<<stretch<<endl;
    cout << "Asym  "<<asymmetry1 << " "<<asymmetry2;
    cout << "  asyE1 "<<asyE1<<"  asyE2 "<<asyE2<<"    Asy = "<<asymmetry<<endl;   
    cout << "energies "<<e1<<"  "<<e2<<"   "<<energy<<endl;
    
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

  Int_t idx;

  // Round off the angle a bit to find the index. 
  
 
  vector<Float_t>::const_iterator iterAngle = upper_bound(angle_row.begin(), angle_row.end(), angle); //search for the first value of angle which is larger than actual

  idx = iterAngle - angle_row.begin()-1;

  if (idx > angle_row.size()-2) idx = -1;

  return idx;
  
}

Int_t hamcPhyPREX::FindEnergyIndex(Float_t energy){
  
  // Round off the energy a bit to find the index.  (E is fixed later by interpolation)
  
  Float_t rounded_energy = round(0.98*energy*20)/20;

  vector<Float_t>::const_iterator iterEnergy  = lower_bound(energy_row.begin(), energy_row.end(), rounded_energy); //search for the first value of energy which is larger or equal than actual 

  UInt_t isize = energy_row.size();
  
  if (rounded_energy > energy_row[isize-1]) return -1;

  return iterEnergy - energy_row.begin(); 
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

  Int_t ldebug=0;

  Float_t pi = 3.1415926; 

  if (nuc != 0 && nuc != 1) {
    cout << "hamcPhyPREX::ERROR: invalid nucleus choice "<<endl;
    return 0;
  }

  if (qsq == 0) {  // wasn't calculated yet.  (qsq is class member)
    qsq = 2*energy*energy*(1-TMath::Cos(angle*pi/180.));
  }

  Float_t calcrsec, mott, form_factor;

  /*Mott cross section for point-like scattering = 
    (alpha*Z* (hc/2pi)*cos(theta/2))^2 / 400E^2(sin(theta/2))^4 */
  
  Float_t halfangle_rad = (angle/2)*(pi/180);
  
  Float_t sin4 = pow(sin(halfangle_rad),4);
  Float_t znuc;
  if (nuc == 0) {
    znuc = 82;
  } else {
    znuc = 6;
  }

  Float_t ascreen = 5e5;  // units: fm (1e-15 m)

  Float_t xscreen = 38.0 / (ascreen * ascreen);

  mott = pow(((znuc*0.197*cos(halfangle_rad))/137),2)/ (xscreen + 400*pow(energy,2)*sin4);

  // qsq comes from hamcKine 

  if (nuc == 0) {  // lead

    vector<Float_t>::const_iterator iterQsq = upper_bound(qsq_row.begin(), qsq_row.end(), qsq); //search for the first value of qsq which is larger or equal than actual   
    int indxQsq = iterQsq - qsq_row.begin(); 

    if (indxQsq == 0) indxQsq = 1; // use the min.

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

  if (ldebug) cout << "CalcCross:  "<<angle<<"  "<<energy<<"   "<<form_factor<<"   "<<
               mott<<"   "<<calcrsec<<endl;

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
  cout << "crsc_table size "<<crsc_table.size()<<"    asym tab size "<<asymmetry_table.size()<<endl;
  
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


void hamcPhyPREX::PrintAsymTable() {

  cout << "\n\n === Asymmetry Table ==== "<<endl;

  for (int str = 0; str < 2; str++) {

    for (int idxE = 0; idxE < 14; idxE++) {

      for (int idxA = 0; idxA < 66; idxA++) {

	cout << "a["<<str<<"]["<<idxE<<"]["<<idxA<<"] = "<< asymmetry_tables[str][idxE][idxA] <<endl;

      }
    }
  }
}
