//  hamcPhyPREX   -- class for the physics of PREX
//  D. Jaunzeikare, R. Michaels  May 2008

// Status Nov 7, 2008
// This is what I got from the summer student Diana.
// Not bad for an 18 y.o. !  
// I still need to run some checks and clean it up,
// but it looks pretty good so far.

#include "hamcPhyPREX.h"
#include "hamcExpt.h"
#include "hamcTarget.h"
#include "hamcEvent.h"
#include "hamcBeam.h"
#include "hamcTrackOut.h"
#include "hamcTrack.h"
#include "hamcInout.h"
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
}


hamcPhyPREX::~hamcPhyPREX() { }


Int_t hamcPhyPREX::Init(hamcExpt* expt) {

  LoadFiles();  // Load the lookup files
  didinit = kTRUE;

  expt->inout->AddToNtuple("crsec",&crsec);
  expt->inout->AddToNtuple("asy",&asymmetry);

  tlen = expt->target->GetRadLength();

  tgtM = expt->target->GetMass();

  /* Initializing Radiative Corrections*/
  ncell1 = 19065; //initializing constants (number of cells)
  ncell2 = 9768;  //@todo This number shouldn't be constant, but stop when 20% are reached

  InitInternalRadCor();
  InitExternalRadCor();

  return 1;
}

Int_t hamcPhyPREX::Radiate(){

  dE_IntBrehm = GenerateDeltaE_internal(); //update hamcPhysics variable. This is happening at each event. 
  dE_ExtBrehm = GenerateDeltaE_external();

  return 1;
}

Int_t hamcPhyPREX::Generate(hamcExpt *expt) {

  Float_t energy = expt->event->beam->GetEnergy();
  Float_t theta = expt->event->trackout[0]->GetTheta();

// Compute the cross section for PREX
  CrossSection(energy, theta, 0);

// Also must compute asymmetry 
   Asymmetry(energy, theta,0);

   return OK;
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
  
  if (debug) cout << "crsec indices "<<energy<<"  "<<angle<<"  "<<indxAngle<<" "<<indxEnergy<<endl;

  Float_t crsc1, crsc2;

  //find the angles above and below the actual angle; used for interpolation
  Float_t angle_upper =angle_row[indxAngle];
  Float_t angle_lower = angle_row[indxAngle-1];

  crsc1 = crsc_tables[stretch][indxEnergy][indxAngle+1]; /*get the cross section value for angle value larger than the actual*/
  crsc2 = crsc_tables[stretch][indxEnergy][indxAngle+2]; //get the cross section value for energy value one below the actual
    // }
  crsec = Interpolate(angle_lower, angle_upper, angle, crsc1, crsc2)/1000;  

  if (debug) {
     cout << "\n Cross section : "<<endl;
     cout << "energy " << energy << " GeV " << "angle " << angle << endl;
     cout <<"Interpolation of crsec"<<crsc1/1000 <<" "<<crsc2/1000<<" "<<crsec<<endl;
  }

  Float_t calccrsec = CalculateCrossSection(energy, angle);
  if (debug) cout <<"calculated cross section " << calccrsec <<endl;

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

  if (indxAngle <= 0) return -1;

  Float_t asymmetry1, asymmetry2;
  Float_t angle_upper =angle_row[indxAngle] ;
  Float_t angle_lower = angle_row[indxAngle-1];
  
  asymmetry1 = asymmetry_tables[stretch][indxEnergy][indxAngle+1];
  asymmetry2 = asymmetry_tables[stretch][indxEnergy][indxAngle+2];
  asymmetry = Interpolate(angle_lower, angle_upper, angle, asymmetry1, asymmetry2);

  int debug=0;

  if (debug) {
    cout<<"angle"<<angle<<endl;
    cout << "Asymmetries"<<asymmetry1 << " "<<asymmetry2 <<" " <<asymmetry<<endl;   }

  return 1;
}


Int_t hamcPhyPREX::FindAngleIndex(Float_t angle){

  vector<Float_t>::const_iterator iterAngle = upper_bound(angle_row.begin(), angle_row.end(), angle); //search for the first value of angle which is larger than actual

  return iterAngle - angle_row.begin();
  
}

Int_t hamcPhyPREX::FindEnergyIndex(Float_t energy){
  
  Float_t rounded_energy = round(energy*20)/20;
  
  vector<Float_t>::const_iterator iterEnergy  = lower_bound(energy_row.begin(), energy_row.end(), rounded_energy); //search for the first value of energy which is larger or equal than actual 

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
  /* Loading Horowitch tables in memory. Data is saved in 3D vector( aVector[stretch][energy][angle]*/

  int stretch = 0;

  vector<vector<Float_t> > crsc_table_temp;
  vector<vector<Float_t> > asymmetry_table_temp;

  LoadHorowitzTable(crsc_table_temp, asymmetry_table_temp,stretch); 
  crsc_tables.push_back(crsc_table_temp);
  asymmetry_tables.push_back(asymmetry_table_temp);

  crsc_table_temp.clear();
  asymmetry_table_temp.clear();

  stretch = 1; 
  LoadHorowitzTable(crsc_table_temp, asymmetry_table_temp,stretch);
  crsc_tables.push_back(crsc_table_temp);
  asymmetry_tables.push_back(asymmetry_table_temp);

  LoadFormFactorTable();
  
  cout<< "hamcPhyPREX:  Tables loaded" <<endl;

  return 1;
}


Float_t hamcPhyPREX::CalculateCrossSection(Float_t energy, Float_t angle) {
 
  Float_t calcrsec, mott, form_factor;

  /*Mott cross section for point-like scattering = 
    (alpha*Z* (hc/2pi)*cos(theta/2))^2 / 400E^2(sin(theta/2))^4 */
  
  Float_t pi = 3.1415926; //@todo Expt class also defines pi
  Float_t halfangle_rad = (angle/2)*(pi/180);
  
  Float_t sin4 = pow(sin(halfangle_rad),4);
  mott = pow(((82*0.197*cos(halfangle_rad))/137),2)/ (400*pow(energy,2)*sin4);
   
  Float_t qsq = CalculateQsq(energy, angle); //@todo move to class hamcKine

  vector<Float_t>::const_iterator iterQsq = upper_bound(qsq_row.begin(), qsq_row.end(), qsq); //search for the first value of qsq which is larger or equal than actual   
  int indxQsq = iterQsq - qsq_row.begin(); 

  Float_t qsq1, qsq2, form_factor1, form_factor2;
  qsq1 =  qsq_row.at(indxQsq);
  qsq2 = qsq_row.at(indxQsq-1);

  form_factor1 = ffsq_row.at(indxQsq);
  form_factor2 = ffsq_row.at(indxQsq-1);
  form_factor = Interpolate(qsq1, qsq2, qsq, form_factor1, form_factor2);

  calcrsec = mott*form_factor; //result is in barn/seradians multiply by 1000 to compare with the value from Horowitchs table where it is milibars/stereadians 
  
  return calcrsec;
}

Float_t hamcPhyPREX::CalculateAsymmetry(Float_t energy, Float_t angle){
  /*A(LR) = (G * Q^2 )/ (4pi*alpha*Sqrt(2))* [1-4sin(theta w)^2 - NeutronFF/ProtonFF ]
  GF = 1.6637 * 10^-5 GeV^-2 Fermi constant
  alpha = 1/137
  sin(theta w)^2 = 0.227 Weinberg Angle

  4*sin^2(theta_W) - N/Z = -1.44
  G / (4 pi alpha sqrt(2)) = 89.4 ppm / GeV^2  
  Asymmetry = - 89.4 * -1.44 * Q^2 = +128.8 Q^2 ppm
  
  */
  //@todo get N and Z from hamcGlobals? or hamcKine?
  Float_t  qsq = CalculateQsq(energy, angle); 
  
  return 128.8*qsq*pow(10.0,-6);
}


Float_t hamcPhyPREX::CalculateQsq(Float_t energy, Float_t angle){
//@todo make this part of hamcKine
  Float_t pi = 3.1415926; 
  Float_t Ebeam = energy;
  Float_t theta_rad = angle*pi/180;
  Float_t  Eprime = Ebeam / (1 + (Ebeam/tgtM)*(1-TMath::Cos(theta_rad)) );
  Float_t qsq = 2*Ebeam*Eprime*(1-TMath::Cos(theta_rad));
  return qsq;
}

Int_t hamcPhyPREX::LoadFormFactorTable(){

  /*Load Lead form factor
   first number is q, fifth is form factor squared 
   vector<Float_t> qsq_row, ffsq_row;  
   q = 0.197 * q;
   qsq = q*q;
  */
  
  FILE *fd;
  char strin[200]; //@todo how large should be? 
  char* filename = "mefcal.pb208_1_25deg_fine.out"; //@todo move to header

  float q, ffsq, qsq;
  float ignore0, ignore1, ignore2, ignore3, ignore4, ignore5, ignore6, ignore7;
  
  fd = fopen(filename, "r");
  if (fd==NULL) {
    printf("ERROR: file %s does not exist \n", filename);
    printf("Bye Bye. \n");
    exit(0);
  }

  while(fgets(strin,1000,fd)!=NULL) {//@todo how large should be 1000? 
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

/*This function reads Horowitchs tables into two 2D vectors, one for crossection, another for asymmetry*/

  vector<float> crsc_row;  //temporary variables
  vector<float> asymmetry_row;

  energy_row.clear();
  asymmetry_row.clear();

  FILE *fd;
  float energy, angle, cross_section, asymmetry;
  float ignore1, ignore2, ignore3;
  char strin[100];

  /*Check wich filename*/
  char* filename;
  if (stretch==0) {
    filename="horpb.dat";
  }  else {
    filename="horpb1.dat";
    }

  fd=fopen(filename, "r");
  if (fd==NULL) {
    printf("ERROR: file %s does not exist \n", filename);
    printf("Bye bye. \n");
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
      sscanf(strin, "%f %f %f %f %f %f",
	     &angle, &cross_section, &ignore1, &asymmetry, &ignore2, &ignore3);
      if (didWriteAngle==0) {
	angle_row.push_back(angle);	
      }
      crsc_row.push_back(cross_section);
      asymmetry_row.push_back(asymmetry);
    }
  }
  crsc_table.push_back(crsc_row);
  asymmetry_table.push_back(asymmetry_row);

  fclose(fd);

  return 1;
}


/*---------------- Radiative Corrections ----------------------*/

//@todo put this in hamcRad
Float_t CalculateIe(Float_t E, Float_t E0, Float_t t){
  /*type=0 -> internal Brehmhstralung                                        
   type=1 -> external Bremhstralung */
   
  Float_t b,psi;
  Float_t x1,x2,x3;

  Float_t z = 82;   // lead  
  Float_t alpha=(1./137.);
  Float_t pi=3.1415926;

  x1 = TMath::Exp((-2./3.)*TMath::Log(z));  // This is Z^(-2/3)                   
  x2 = TMath::Exp((-1./3.)*TMath::Log(z));
  psi = TMath::Log(1440*x1)/TMath::Log(183*x2);
  b = (4./3.)*(1 + (1./9.)*(( (z+1)/(z+psi) ) / (TMath::Log(183*x2))));

  Float_t qq,qsq,msq;
  Float_t qf=0.45; // inverse fermis (Q)                                         
  Float_t me=0.511;  // mass electron (MeV)                                      
  qq = qf*0.197;
  qsq = qq*qq;
  msq = (me/1000)*(me/1000);

  // cout << "qsq, msq, teq "<<qsq<<"  "<<msq<<"  "<<t<<endl;
  /*t cant be moved outside of this function, because alpha b are defined here*/ 
  //  cout<<"tlen" <<t<<endl;
 if(t==0){
    t = (alpha/(b*pi)) * (TMath::Log(qsq/msq) - 1); //@todo shoud qsq be taken from hamcExpt? 
 }
 // cout <<"tlen new "<<t<<endl;
  x1 = b*t/(E0-E);
  x2 = E/E0 + (3./4.)*(((E-E0)/E0)*((E-E0)/E0));
  x3 = TMath::Exp(b*t*TMath::Log(TMath::Log(E0/E)));

  return x1*x2*x3;
}


Int_t hamcPhyPREX::InitRadCor(Float_t *elist, Int_t ncell1, Float_t t, TH1F *hb){
 

  Int_t Npt= 2000;
  Float_t E0 = 1.0; //@todo take it from hamcExpt; 

  //TH1F *hhh = new TH1F("h2b","External Bremhstrahlung Ncell",Npt/10,-0.1,1.04*E0);

  Int_t lout=0;
  Float_t Ie, Prob, Ptot1, xnorm1;
 
  Float_t yfact=1e7;
  Float_t E, dE;
 
  Int_t nybin=1000;

  Float_t ycell = 0.05*yfact/((Float_t)nybin); // 0.05*1e7/1000 = 500           
  
  Float_t ncell;

  Int_t idx1;

  // Float_t y1hi=0.1*yfact; //0.1*1e7  //for histograms


  dE = E0/((Float_t)Npt);  // interval                                            
  E=E0;  // Initialize                                                            
  Ptot1=0; // Total probability?
  xnorm1 = 0.612015; //?                                                        
 
  Int_t nct1=0;
  idx1=0; 

  for (Int_t i=0; i<Npt; i++) {
    E = E-dE;
    if (E<0) continue;

    //For each energy calculate function value divide it by yfact to find out how many cells fit vertically and then for each cell add to elist that energy  
    // cout<<"t to pass to CalculateIe "<<t<<endl;
    Ie=CalculateIe(E, E0, t);//calculate Ie for internal Bremshtralung         
    Prob = Ie*dE / xnorm1; //=Ie*dE/ 0.612015 ?Probability                        
    if (i==0) cout << "(straggling) Prob0 "<<yfact*Prob<<endl;
    Ptot1 += Prob; //total probability                                          
    Prob = yfact*Prob; //1e7*Prob                                                
    ncell = Prob/ycell; //Nth cell equals Prob/500                               
    nct1 += ((Int_t)ncell);//ncell count?         

    if (nct1>=ncell1) {
      //      cout << "ERROR: ncell1 too small !!"<<endl;
    } else {
      Int_t nlim=nct1-idx1+1;
      for (Int_t jj = 0; jj<nlim; jj++) {
	elist[idx1] = E;
	//cout << " E "<<E;
	idx1++;
      }
    }

    if(lout==1) cout << "Prob "<<E<<"  "<<Prob<<"   "<<ncell<<endl;

    // h1->Fill(E,Prob);                                                        
    // h1a->Fill(E,Prob);                                                       
    hb->Fill(E,ncell);//@todo ?
    }
  return 1;
  //hhh->Write("hhh");
}

Int_t hamcPhyPREX::InitInternalRadCor(){
  cout<<"Starting InitInternalRadCor()"<<endl;
    Float_t Npt = 2000;
    Float_t E0 = 1.0;
   elistInternal= new Float_t[ncell1];  
   // string TitleString = new string("Internal Bremsstrahlung NCell, tlen = ");
   TH1F *h1b = new TH1F("h1b","Internal Bremsstrahlung NCell, t=0.1",(Int_t)(Npt/10),-0.1,1.04*E0);
   InitRadCor(elistInternal, ncell1, tlen, h1b); 
   //   InitRadCor(elistInternal, ncell1, 0.1, h1b); @todo if tlen doesn't seem to fetch the right value, set here manually.
   //   The current directory (root) is not associated with a file. The object (h1b) has not been written.
 
  h1b->Write("h1b");
   return 1;
}

Int_t hamcPhyPREX::InitExternalRadCor(){
  cout<<"Starting InitExternalRadcor()"<<endl;
  Float_t Npt = 2000;
    Float_t E0 = 1.0;
  elistExternal = new Float_t[ncell2];
  TH1F *h2b = new TH1F("h2b","External Bremsstrahlung Ncell",(Int_t)(Npt/10),-0.1,1.04*E0);
  InitRadCor(elistExternal, ncell2,0, h2b); //zero means to use teqviv.
  h2b->Write("h2b"); // The object (h2b) has not been written.

  return 1;
}

Float_t hamcPhyPREX::GenerateDeltaE_internal(){
 
  Int_t idx = (Int_t)(ncell1*gRandom->Rndm(1));
  return elistInternal[idx];

}

Float_t hamcPhyPREX::GenerateDeltaE_external(){

  Int_t idx = (Int_t)(ncell2*gRandom->Rndm(1));
  return elistExternal[idx];

}

/*Int_t hamcPhyPREX::TestRadiate(){

  Int_t Npt = 2000;
  Float_t E0 = 1.0;

  TH1F *hbrem1 = new TH1F("hbrem1","Straggling",Npt/2,-0.1,1.04*E0);

  for (Int_t i = 0; i<50000; i++) {

    int idx1 = (Int_t)(ncell1*gRandom->Rndm(1));

    Float_t num = elistInternal[idx1];
    // cout <<"num "<<num<<endl;
    if (idx1<0 || idx1>ncell1) {
      cout << "ERROR:  idx1 out of range "<<endl;
    } else {
      //hbrem.Fill(elist1[idx1]);                                                 
      hbrem1->Fill(num);
    }
  }
  hbrem1->Write("hbrem1");
  
  return 1;
}
*/
