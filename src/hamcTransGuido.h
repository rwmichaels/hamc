#ifndef ROOT_hamcTransGuido
#define ROOT_hamcTransGuido

//  hamcTransGuido -- Guido's functions to transport to focal plane
//  R. Michaels  Mar 2009

#include "Rtypes.h"
#include "hamcTrans.h"
#include <vector>
#include <iostream>
#include <string>
#include <map>

class hamcSpecHRS;

class hamcTransGuido : public hamcTrans {

  public:

     hamcTransGuido();
     virtual ~hamcTransGuido(); 
     Int_t Init(hamcSpecHRS *spect);
     void Print();

// Transform 'trk' to 'where'.
     Int_t TransForm(hamcTrack *trk, Int_t where) const;  // modifies trk

// Transform a transport vector
     Int_t TransForm(hamcTransVect *tvect, Int_t where) const; // modifies tvect

  private: 

     hamcTransGuido(const hamcTransGuido& trans);
     hamcTransGuido& operator=(const hamcTransGuido& trans);

// Parameters of the transport
     //Float_t xDth,xDthth,xDph,xDphph,yDth,yDph;
     //Float_t yDphph,thDx,phDx,phDxx,phDy;
     Float_t  xDph,xDdpp,xDdppdpp,xDdppdppdppdpp,xDth,xDththth,xDthththth,xDththththth,xDthdpp,xDthththph,xDththththththph,xDthdppdpp,xDthththdppdpp,xDthththdppdppdpp,xDththththdppdppdpp,xDphdpp,xDphdppdppdpp,xDthphdpp,xDthy0y0,xDthy0y0y0,xDdppy0,xDy0;
     Float_t yDph,yDdpp,yDth,yDthth,yDththth,yDphdpp,yDphphdpp,yDy0;
     
#ifndef NODICT
ClassDef (hamcTransGuido, 0)   // Guido's Model for HRS Transport
#endif
};

#endif



   
