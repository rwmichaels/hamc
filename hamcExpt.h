#ifndef ROOT_hamcExpt
#define ROOT_hamcExpt

//  hamcExpt   -- abstract base class for experiment
//  R. Michaels  June 2008

#include "Rtypes.h"
#include "hamcSpecHRS.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>

#define IELECTRON 10
#define IPROTON   11
#define IPION     12
#define melectron 0.000511
#define mpion     0.140
#define mproton   0.938
#define PI        3.1415926

#define OK     1
#define ERROR -1

class hamcTarget;
class hamcPhysics;
class hamcEvent;
class hamcInout;

class hamcExpt {

  public:

     hamcExpt(std::string sname);
     virtual ~hamcExpt()=0;    // abstract class

     virtual Int_t Init();
     virtual Int_t Init(std::string sfile);
     virtual Int_t InitInput(std::string sfile);
     virtual Int_t Run(Int_t maxev=1000000);
     virtual Int_t Summary();
     Int_t GetNumSpectrom() const { return (Int_t)spectrom.size(); };
     hamcSpecHRS* GetSpectrom(Int_t i) const;
 // Decided to make these public instead of using "get" functions
     hamcTarget  *target;
     hamcPhysics *physics;
     hamcEvent   *event;
     hamcInout   *inout;
     Int_t iteration;    // which iteration we are on.
     std::string name;

     std::string GetSetupFile() { return setupfile; };

  protected:

     virtual Int_t ReadSetup();
     virtual Int_t SetIterate();
     virtual void  Analysis()=0;  // must define in the inheriting class

     Int_t ChkSpIndex(Int_t i) const;

     std::vector<hamcSpecHRS* > spectrom;

     Int_t  numevent, numiter;
     std::string setupfile;
     Bool_t didinit;
     std::vector<std::string> setupstr;

     hamcExpt(const hamcExpt& expt);
     hamcExpt& operator=(const hamcExpt& expt);
     
  private: 

#ifndef NODICT
ClassDef (hamcExpt, 0)   // An experiment
#endif


};

#endif



   
