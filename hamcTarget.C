//  hamcTarget   -- Target
//  R. Michaels  June 2008

#include "hamcTarget.h"
#include "hamcExpt.h"
#include "Rtypes.h"
#include "TRandom.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;


#ifndef NODICT
ClassImp(hamcTarget)
#endif

hamcTarget::hamcTarget(string tgt_name) : name(tgt_name),did_init(kFALSE)
{
  overall_length=0;
  effective_length=0;
  radiation_length=0;
  width=0;
  weighted_anum=0;
  weighted_znum=0;
  weighted_density=0;
  mass=0.938;   // a default
  zscatt=0;
  material_index=0;
}

hamcTarget::~hamcTarget() {
}

Int_t hamcTarget::Setup() {

  if (did_init) return OK;

  if (components.size() == 0) {
    cout << "hamcTarget::Setup: WARNING: No target components !"<<endl;
  } else {

    Float_t weisum=0;

    for (Int_t i=0; i<(Int_t)components.size(); i++) {

       hamcTgtSlab *slab = components[i];
       overall_length   += slab->GetLen(); 
       radiation_length += slab->GetLen()/slab->GetRadLen();
       Float_t weight = slab->GetDensity()*slab->GetLen();
       weisum           += weight;
       effective_length += weight*slab->GetEffLen();
       weighted_anum    += weight*slab->GetA();
       weighted_znum    += weight*slab->GetZ();
       weighted_density += weight*slab->GetDensity();
    }

    if (weisum == 0) {
       cout << "hamcTarget::ERROR: zero weights !"<<endl;
    } else {
       effective_length = effective_length/weisum;
        weighted_anum   = weighted_anum/weisum;
       weighted_znum    = weighted_znum/weisum;
       weighted_density = weighted_density/weisum;
    }

    Float_t zpos;
    zpos = -1 * overall_length/2;  // z=0 is center of overall target

    for (Int_t i=0; i<(Int_t)components.size(); i++) {

       zpos += (components[i]->GetLen())/2;  // 1st half

       components[i]->PutZloc(zpos);

       zpos += (components[i]->GetLen())/2;  // 2nd half


    }

    did_init = kTRUE;

  }

  FindMtlIndex(0);

  return OK;

}

Int_t hamcTarget::Zscatt() {

   zscatt = overall_length * (-0.5 + gRandom->Rndm());

   return FindMtlIndex(zscatt);

}			      


Int_t hamcTarget::FindMtlIndex(Float_t zloc) {
// Find the index of the material for this Z location.
// The index tracks with (A,Z) and is use to lookup cross-section.
// The 'mass' is also looked up, as it is needed for kinematics.

    Float_t zdiff;

    material_index = -1;
    mass = 0; 

    for (Int_t i=0; i<(Int_t)components.size(); i++) {

       hamcTgtSlab *slab = components[i];
       zdiff = zloc - slab->GetZloc();
       if (zdiff < 0) zdiff = -1.*zdiff;
// within 1/2 of length, the 0.502 is to account for roundoff
       if (zdiff < 0.502*slab->GetLen()) {
         material_index = slab->GetIndex();
         mass = slab->GetMass();
         return OK;
       }
    }

// This would be a pretty serious error.  Should fix it if it happens
    cout << "hamcTarget::FindMtlIndex:WARNING: No component found at this Z"<<endl;
 
    return OK;

}


void hamcTarget::Print() {

  cout << "hamcTarget::Print"<<endl;
  cout << "Target name : "<<GetName()<<endl;
  cout << "Length "<<GetLength()<<"  Eff. Length "<<GetEffLength()<<endl;
  cout << "Radiation Length = "<<GetRadLength()<<endl;
  cout << "<A> = "<<GetA()<<"   <Z> = "<<GetZ()<<endl;
  cout << "<density> = "<<GetDensity()<<endl;
  cout << "present Z scatter point "<<GetZScatt()<<endl;
  cout << "mass of Z scatter point = "<<GetMass()<<endl;
  cout << "Number of materials "<<GetNumMtl()<<endl;
  for (Int_t i=0; i<GetNumMtl(); i++) {
    cout << "Material "<<i<<"  A = "<<GetMtlA(i)<<"  Z = "<<GetMtlZ(i)<<endl;
    cout << "Z location "<<GetMtlZloc(i);
    cout <<"    RL = "<<GetMtlRadLen(i)<<"   rho = "<<GetMtlDensity(i)<<endl;
  }

}


Float_t hamcTarget::GetMtlZloc(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetZloc();

}


Float_t hamcTarget::GetMtlRadLen(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetRadLen();

}

Float_t hamcTarget::GetMtlDensity(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetDensity();

}

Float_t hamcTarget::GetMtlA(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetA();

}

Float_t hamcTarget::GetMtlZ(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetZ();

}


Int_t hamcTarget::CheckIndex(Int_t index) {

   if ( !did_init ) {
     cout << "hamcTarget::CheckIndex:ERROR: hamcTarget not initialized"<<endl;
     return ERROR;
   }
 
   if (index < 0 || index > ((Int_t)components.size())-1) {
     cout << "hamcTarget::CheckIndex:ERROR: index out of range !"<<endl;
     return ERROR;
   }

   return OK;

}





