
//  hamcSpecHRS   -- HRS Spectrometer
//  R. Michaels  June 2008

#include "hamcSpecHRS.h"
#include "hamcInout.h"
#include "hamcExpt.h"
#include "hamcTransMat.h"
#include "hamcTHRSTrans.h"
#include "hamcTransLerHRS.h"
#include "hamcTransLerColdSeptum.h"
#include "hamcTransLerWarmSeptum.h"
#include "hamcTransLer4deg.h"
#include "hamcTransGuido.h"
#include "hamcTrans.h"
#include "hamcAperture.h"
//#include "hamcDetector.h"
#include "Rtypes.h"
#include <string>
#include <vector>

using namespace std;

#ifndef NODICT
ClassImp(hamcSpecHRS)
#endif

hamcSpecHRS::hamcSpecHRS(Int_t which, Float_t pmom, Float_t angle) : P0(pmom), central_angle(angle) {
  which_spectrom = which;  // LEFT, RIGHT 
  name = "HRS";
  desc = "Hall A High Resolution Spectrometer";
  P0_sigma = P0*1e-4;       // default resolution
  collim_distance = 1.1825;   // <- meters, used to be 1.24
  // Formula to get quad values:  Take HRS current (A) and apply this
  // factor:  1.458e-3 for Q1 and 2.025e-3 for Q2,Q3
  // Signs are + - -
  // Ryan's defaults (Nov, 2017)
  quad1 = 0.096255;    //  490 A
  quad2 = -0.131739;   // 429 A
  quad3 = -0.170480;   // 579 A
  // Here is for Tune B (PREX-I) with a focus 1.43m after VDC1.
  //  quad1 = 0.116;   // 582 A
  // quad2 = -0.128;    // 462 A
  // quad3 = -0.158;    // 570 A
  // And I think these dipole parameters never change
  dipk1 = 0.050178;
  dipk2 = 0.037056;
  fTune = THRSTrans::kPREX;

  // set collim_distance (numbers from Kiad's pointing diagrams), 13Dec10, rupesh
  // collim_distance below are for hapIII
  if (which_spectrom == LEFTHRS) collim_distance = 1.1825;   // <- meters LHRS
  else if(which_spectrom == RIGHTHRS) collim_distance = 1.1759;   // <- meters RHRS

  sept_choice   = noseptum; 
  trans_choice  = tmatrix;  
  collim_choice = nocollim;
  use_guido = kFALSE;
  transport = 0;
  tguido = 0;
}

// everything below is from Diancheng.
// acceptance cuts have been fine tuned based on 
// VDC x,y,theta,phi,dp  
// 13Dec10, rupesh

hamcSpecHRS::~hamcSpecHRS() {
  if (transport) delete transport;
  for (vector<hamcSpecBrk*>::iterator ith = break_point.begin();
      ith != break_point.end(); ith++) delete *ith;
}

void hamcSpecHRS::UseCollimator(void) {
  collim_choice = reg_coll;
  cout << "Using regular collimator "<<endl;
}

void hamcSpecHRS::UsePaulColl(void) {
  collim_choice = paul_coll;
}

void hamcSpecHRS::UseAngleColl(void) {
  // Empirical angle cut.
  collim_choice = angle_coll;
}

void hamcSpecHRS::UseHRSOnly() {
  sept_choice = noseptum;
}

void hamcSpecHRS::Use4degSeptum() {
  sept_choice = sept4deg;
}

void hamcSpecHRS::UseWarmSeptum() {
  sept_choice = warmsept;
}

void hamcSpecHRS::UseColdSeptum() {
  sept_choice = coldsept;
}

void hamcSpecHRS::UseTHRSTrans() {
  trans_choice = hrstrans;
}

void hamcSpecHRS::UseMatrixTrans() {
  trans_choice = tmatrix;
}

void hamcSpecHRS::UseLeroseTrans() {
  trans_choice = tlerose;
}

void hamcSpecHRS::UseGuidoTrans() {
  use_guido = kTRUE;
}


Int_t hamcSpecHRS::Init(hamcExpt *expt) {

// Using hamcInout and string parser to obtain setup info.

// The user class (e.g. hamcExptPREX) may define the defaults,
// but they can also be defined (over-ride) here with the setup file.

   cout << "Init hamcSpecHRS : "<<endl;

   hamcStrParser parser;

   parser.Load(expt->inout->GetStrVect("hrs_setup"));
   //   parser.Print();
   if (parser.IsFound("noseptum")) {
     UseHRSOnly();
   }   
   if (parser.IsFound("thrstrans") || parser.IsFound("hrstrans")) {
     cout << "hamcSpecHRS: Using THRSTrans 2nd order matrix model"<<endl;
     UseTHRSTrans();
   } else {
     cout << "Warning:  Should be using thrstrans, the fortran transport"<<endl;
     cout << "  becoming obsolete ! "<<endl;
   }  
   if (parser.IsFound("coldseptum")) {
     cout << "hamcSpecHRS: Using cold septum"<<endl;
     UseColdSeptum();
   }   
   if (parser.IsFound("sept4deg")) {
     cout << "hamcSpecHRS: Using the 4-degree septum"<<endl;
     Use4degSeptum();
   }   
   if (parser.IsFound("warmseptum")) {
     cout << "hamcSpecHRS: Using warm septum"<<endl;
     UseWarmSeptum();
   }   
   if (parser.IsFound("usecollimator")) {
     cout << "hamcSpecHRS: Using collimator"<<endl;
     UseCollimator();
   }   
   if (parser.IsFound("usepaulcollim")) {
     cout << "hamcSpecHRS: Using Paul's composite collimator"<<endl;
     UsePaulColl();
   }   
   if (parser.IsFound("useempiricalangle")) {
     cout << "hamcSpecHRS: Using empirical angle collimation"<<endl;
     UseAngleColl();
   }   
   if (parser.IsFound("usematrix")) {
     cout << "hamcSpecHRS: Using matrix transport"<<endl;
     UseMatrixTrans();
   }   
   if (parser.IsFound("uselerose")) {
     cout << "hamcSpecHRS: Using LeRose functions"<<endl;
     UseLeroseTrans();
   }   
   if (parser.IsFound("useguido")) {
     // Guido's parameterziation has no acceptance criteria.
     // Therefore we use this as an alternative to transport
     // but must use another model for acceptance.
     cout << "hamcSpecHRS: Using Guido's parameterization"<<endl;
     UseGuidoTrans();
     tguido = new hamcTransGuido();
     tguido->Init(this);
   }   

   parser.Load(expt->inout->GetStrVect("quad_fields"));
   // Field gradients (Tesla/m)
   if (parser.IsFound("q1")) {
     quad1 = parser.GetData(); 
   }      
   if (parser.IsFound("q2")) {
     quad2 = parser.GetData(); 
   }      
   if (parser.IsFound("q3")) {
     quad3 = parser.GetData(); 
   }      
   if (parser.IsFound("tune")) {
     quad3 = parser.GetData(); 
   }      
   parser.Load(expt->inout->GetStrVect("quad_tune"));
   if (parser.IsFound("std")) fTune = THRSTrans::kStd;
   if (parser.IsFound("prex")) fTune = THRSTrans::kPREX;
   if (parser.IsFound("crex")) fTune = THRSTrans::kCREX;

   collim2_radlen1 = 0;

   collim2_a = 0;
   collim2_z = 0;
   collim2_t = 0;

   parser.Load(expt->inout->GetStrVect("spreader_collim"));
   if (parser.IsFound("radlen1")) {
     collim2_radlen1 = parser.GetData(); 
   }      
   if (parser.IsFound("a")) {
     collim2_a = parser.GetData(); 
   }      
   if (parser.IsFound("z")) {
     collim2_z = parser.GetData(); 
   }      
   if (parser.IsFound("t")) {
     collim2_t = parser.GetData(); 
   }      



   BuildSpectrom();

// For test mode we may turn off acceptance cuts in some places.
   vector<string> sdata; 
   sdata = expt->inout->GetStrVect("HRS_acceptoff");
   if (sdata.size()>=1) {
     cout << "WARNING: turning off acceptance cuts in some places "<<endl;
     for (Int_t isd = 0; isd < (Int_t)sdata.size(); isd++) {
       if (isd < (Int_t)break_point.size()) 
         break_point[isd]->DefineCut(atoi(sdata[isd].c_str()));
     }
   }

   if (transport) transport->Init(this);
   
   return OK;

}

Int_t hamcSpecHRS::BuildSpectrom() {

// Build the HRS transport and acceptance model
// depending on the choices.

  AddBreakPoint(ITARGET);

  if (IsMatrixTrans()) {   // Transport matrix

     transport = new hamcTransMat();

     if (IsCollimated()) {
        if (IsPaulCollim()) {
 	   AddBreakPoint(ICOLLIM2);
	} else {
	   if (IsAngleCollim()) {
   	      AddBreakPoint(ICOLLIM3);
	   } else {
              AddBreakPoint(ICOLLIM);
	   }
        }
     }
     AddBreakPoint(IFOCAL); 

     return OK;
  } 

  if (IsTHRSTrans()) { // Seamus's 2nd order TRANSPORT model

     transport = new hamcTHRSTrans();

     if (IsCollimated()) {
        if (IsPaulCollim()) {
 	   AddBreakPoint(ICOLLIM2);
	} else {
	   if (IsAngleCollim()) {
   	      AddBreakPoint(ICOLLIM3);
	   } else {
              AddBreakPoint(ICOLLIM);
	   }
        }
     }
     if (IsWarmSeptum() || IsColdSeptum()) {
        AddBreakPoint(ISEPTIN);   // septum input
        AddBreakPoint(ISEPTOUT);  // setpum out
     } 

     AddBreakPoint(IQ1IN);
     AddBreakPoint(IQ1EXIT);
     AddBreakPoint(IQ2IN);
     AddBreakPoint(IQ2EXIT);
     AddBreakPoint(IDIPIN);
     AddBreakPoint(IDIPEXIT);
     AddBreakPoint(IQ3IN);
     AddBreakPoint(IQ3EXIT);
     AddBreakPoint(IFOCAL);

  }

#ifdef USEFORTRAN

  if (IsLeroseTrans()) {   // LeRose functions

    if (Is4degSeptum()) {
      transport = new hamcTransLer4deg();  // 4-degree Septum
    }
    if (IsWarmSeptum()) {
      transport = new hamcTransLerWarmSeptum();  // Warm Septum
    } 
    if (IsColdSeptum()) {
      transport = new hamcTransLerColdSeptum();  // Cold Septum
    } 
    if (!IsWarmSeptum() && !IsColdSeptum() && !Is4degSeptum()) {    // HRS w/o Septum
      transport = new hamcTransLerHRS();
    }


    if (IsCollimated()) {
        if (IsPaulCollim()) {
 	   AddBreakPoint(ICOLLIM2);
	} else {
	   if (IsAngleCollim()) {
   	      AddBreakPoint(ICOLLIM3);
	   } else {
              AddBreakPoint(ICOLLIM);
	   }
        }
     }

    if (IsWarmSeptum() || IsColdSeptum() || Is4degSeptum())  {

       AddBreakPoint(ISEPTIN);   // septum input
       AddBreakPoint(ISEPTOUT);  // setpum out
    }
        
    AddBreakPoint(IQ1EXIT);
    AddBreakPoint(IDIPIN);
    AddBreakPoint(IDIPEXIT);
    AddBreakPoint(IQ3IN);
    AddBreakPoint(IQ3EXIT);
    AddBreakPoint(IFOCAL);
     
    return OK;
  }

#else

  if (IsLeroseTrans()) {   // LeRose functions
    cout << "hamcSpecHRS::WARNING: No Fortran functions allowed in this version"<<endl;
    cout << "Therefore, no LeRose functions for Transport defined"<<endl;
    return -1;
  }
#endif

// If there are Guido functions, add the choice here.

  return OK;

}


void hamcSpecHRS::AddBreakPoint(Int_t where) {

// Add a break point with an aperture to define acceptance
// that depends on "where"
// The numbers here are just a dream at this point. (need to fix).

  Int_t idx;

  printf("*************************   where = %d\n", where );

  switch(where) {

    case ITARGET:
      break_point.push_back(new hamcSpecBrk(where));
      break;

    case ICOLLIM:
      if (IsWarmSeptum() || IsColdSeptum()) {
        if (IsWarmSeptum()) {
  	   break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.125,0.12,-0.12,0.12)));
 	} else {
 	  break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.30,-0.2156,-0.10,0.10)));  // cold septum
 	}
       } else { // HRS alone  collimator: (horiz)62.9 mm x  (vert)121.8 mm 
       //	break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.0609,0.0609,-0.03145,0.03145)));  // collimator entrance
       //	break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.06485,0.06485,-0.0334,0.0334)));  // collimator exit
	
	// blank out Diancheng's collim below
 	//	break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.06485,0.06485,-0.0285,0.0334)));  // collimator exit

 	// tmpoffsets are offsets to collimator, used to move hapIII collimator
 	// correct signs now
 	// 13Oct10, rupesh
 	Float_t tmpoffsetx = 0.0, tmpoffsety=0.0 ;

 	// unblank below for hap III
// 	// 0.67 mm towards beamline for LHRS
// 	if (which_spectrom == LEFHRS) tmpoffsetx = 0.00067;  //0.67mm
// 	// 1.04 mm towards beamline for RHRS
// 	else if (which_spectrom == RIGHTHRS) tmpoffsetx = -0.00104; //1.04mm
	
// 	cout << "The collimator offset in x is "<< tmpoffsetx << endl;

         // tmpsizeoffset is used to change collim size for HRS acceptance tests for hapIII
 	Float_t tmpsizeoffset = 0.00;  // +0.005 changes collim size by - 5mm	

         break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.0609+tmpoffsety+tmpsizeoffset,
  								 0.0609-tmpoffsety-tmpsizeoffset,
     								 -0.03145-tmpoffsetx+tmpsizeoffset
 								 ,0.03145-tmpoffsetx-tmpsizeoffset)));  

       }
       idx = break_point.size()-1;
       break_point[idx]->aperture->SetCenter(0,0);
       break;

// Collim2 is in the Q1 coordinate system (center = 0,0).
     case ICOLLIM2:
       if (IsWarmSeptum()) { 
        break_point.push_back(new hamcSpecBrk(where, new hamcPaulCollim(
          0.032, 0.041,  0.20, 0.229,  // A_T hole (low, right, and R,C of arc)
  	 0.205, 0.145,  0.20, 0.229,     // outer, inner circles
          0.117, 0.04,       // top, right
          0.1474, -1.88)));  // Champhor line.
         idx = break_point.size()-1; 
         break_point[idx]->aperture->DefineMaterial(0,collim2_a, collim2_z, collim2_t);
 	printf("%f %f %f\n", collim2_a, collim2_z, collim2_t);
       }
       break;

// Collim3 is the empirical angle collimation which is a tighter cut than COLLIM2.
     case ICOLLIM3:
       if (IsWarmSeptum()) { 
 	break_point.push_back(new hamcSpecBrk(where, new  hamcAngleCollim()));
         idx = break_point.size()-1; 
       }
       break;

     case ISEPTIN:
       if (IsWarmSeptum()) break_point.push_back(new hamcSpecBrk(where, new hamcBox(0.088,0.382,-0.12,0.12)));
       if (IsColdSeptum()) break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.1486,-0.0887,-0.110,0.110)));
       break;

     case ISEPTOUT:
       if (IsWarmSeptum()) break_point.push_back(new hamcSpecBrk(where, new hamcBox(0.088,0.382,-0.12,0.12)));
       if (IsColdSeptum()) break_point.push_back(new hamcSpecBrk(where, new hamcBox(-0.3485,-0.2156,-0.110,0.110)));
       if (Is4degSeptum()) {
          break_point.push_back(new hamcSpecBrk(where, new hamcBox(0.088,0.382,-0.12,0.12)));
       }
       break;

     case IDIPIN:
       // Needs some logic here for THRSTrans ... the -5.22 to -4.98 cut is for LeRose only
       if (IsLeroseTrans()) {   // LeRose functions
         if (IsWarmSeptum() || IsColdSeptum()) {
           break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-5.22, -4.98, -0.1924, -0.1924)));
           break;
         }
         if (Is4degSeptum()) {
           break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-6.19, -5.9, -0.161, -0.1924)));
           break;
         }
         break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-0.4, 0.4, 0.125, -0.0186)));  // standard HRS
       } else {  // IsTHRSTrans() or other model
         break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-0.4, 0.4, 0.125, -0.0186)));  
       }
       break;

     case IDIPEXIT:
       if (IsWarmSeptum() || IsColdSeptum()) {
         break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-0.462, 0.462, 0.125, -0.0161)));
         break;
       }
       if (Is4degSeptum()) {
         break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-0.462, 0.462, 0.125, -0.0161)));
         break;
       }
       break_point.push_back(new hamcSpecBrk(where, new hamcTrapezoid(-0.4, 0.4, 0.125, -0.0186)));  // standard HRS

       break;

     case IQ1EXIT:
       break_point.push_back(new hamcSpecBrk(where, new hamcCircle(0.1492)));
       break;

     case IQ3IN:
       break_point.push_back(new hamcSpecBrk(where, new hamcCircle(0.3)));
       break;

     case IQ3EXIT:
       break_point.push_back(new hamcSpecBrk(where, new hamcCircle(0.3)));

       break;

     case IFOCAL:
       break_point.push_back(new hamcSpecBrk(where));
       break;

     default:
       cout << "hamcSpecHRS::WARNING: undefined location "<<where<<endl;
  }
}

 hamcAperture* hamcSpecHRS::Aperture(Int_t idx) {
   if (ChkIdx(idx) == ERROR) return 0;
   return break_point[idx]->aperture;
 }

 Int_t hamcSpecHRS::GetNumBrk() {
   return break_point.size();
 }

 Int_t hamcSpecHRS::ChkIdx(Int_t idx) {
   if (idx < 0 || idx > (Int_t)break_point.size()-1) {
     cout << "hamcSpecHRS::Warning: index out of range"<<endl;
     return ERROR;
   }
  return OK;
 }

void hamcSpecHRS::Print() {
  cout << endl << "hamcSpecHRS  print "<<endl;
  cout << "name = "<<name<<endl;
  cout << "description = "<<desc<<endl;
  cout << "P0 = "<<GetP0() << "    P0 sigma "<<GetP0Sigma()<<endl;
  cout << "Angle = "<<GetScattAngle()<<endl;
  cout << "Septum choice "<<sept_choice<<endl;
  cout << "Collimator choice "<<collim_choice<<endl;
  cout << "Collim distance = "<<GetCollimDist()<<endl;
  cout << "Number of break points = "<<break_point.size()<<endl;
  cout << "Quads "<<quad1<<"  "<<quad2<<"  "<<quad3<<"   dipole "<<dipk1<<"  "<<dipk2<<endl;
  cout << "Tune "<<fTune<<endl;
  for (Int_t i=0; i<(Int_t)break_point.size(); i++) {
    break_point[i]->Print();
  }
}
