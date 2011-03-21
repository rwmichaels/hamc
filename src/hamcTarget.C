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

  // Use zoffset to move the target along the central beam line to the dump

  zoffset = 0.0; 

  // below are happex III offset values, rupesh 13Dec10
  // LHRS, need to move the target downstream along the central beam axis
  //  zoffset = 0.004183; // in meters LHRS,  = 1.02/sin(0.2463) mm 
  // RHRS, need to move the target upstream along the central beam axis
  //zoffset = -0.013058; // in meters RHRS, = 3.09/sin(0.2389) mm
}

hamcTarget::~hamcTarget() {
}

Int_t hamcTarget::Setup() {

  if (did_init) return OK;

  if (components.size() == 0) {

    cout << "hamcTarget::Setup: WARNING: No target components !"<<endl;
    exit(0);

  } else {

    Float_t weisum=0;

    for (Int_t i=0; i<(Int_t)components.size(); i++) {

       hamcTgtSlab *slab = components[i];
       overall_length   += slab->GetLen(); 
       Float_t weight = slab->GetDensity()*slab->GetLen();
       weisum           += weight;
       effective_length += weight*slab->GetEffLen();
       weighted_anum    += weight*slab->GetA();
       weighted_znum    += weight*slab->GetZ();
       weighted_density += weight*slab->GetDensity();
    }

    if (weisum == 0) {
       cout << "hamcTarget::ERROR: zero weights !"<<endl;
       exit(0);
    } else {
       effective_length = effective_length/weisum;
       weighted_anum    = weighted_anum/weisum;
       weighted_znum    = weighted_znum/weisum;
       weighted_density = weighted_density/weisum;
    }

    Float_t zpos;

    // -ve is away from the spec.
    // zoffset added to move target rel to the HRS, 14Oct10 
    zpos = -1 * overall_length/2 + zoffset;  // z=0 is center of overall target

    cout << "The target is offset from the hall center by \t" << zoffset << endl;
//     cout << "overall_length\t" << overall_length << endl;
//     cout << "zpos\t" << zpos << endl;

 // A "weighted radiation_length" is nonsense (this is not a mix),
 // so we do the following to get an approximate total radlength.
    zscatt = 0;
    Float_t x1 = GetRadIn(); 
    Float_t x2 = GetRadOut(); 
    radiation_length = x1 + x2;


    cout << "Weighted radiation length  "<<radiation_length<<endl;

    for (Int_t i=0; i<(Int_t)components.size(); i++) {

      zpos += (components[i]->GetLen())/2;  // 1st half

      components[i]->PutZloc(zpos);
      
      zpos += (components[i]->GetLen())/2;  // 2nd half

    }

    did_init = kTRUE;

    Print();

  }

  FindMtlIndex(0);

  return OK;

}

Int_t hamcTarget::Zscatt() {

  // zoffset added to move target rel to the HRS, 14Oct10 
  zscatt = overall_length * (-0.5 + gRandom->Rndm()) + zoffset;
  //  cout << "zoffset\t" << zoffset << endl;
  //  cout << "zscatt\t" << zscatt << endl;

  return FindMtlIndex(zscatt);
}			      

Float_t hamcTarget::GetRadIn() {
// material_index = index of material where scattering occured
// Returns the fractional radiation length going in to scattering point.

  if (material_index == -1) return 0;
  Float_t rlen = 0;
  Float_t zloff = 0;
  Float_t zdist;
  
  for (Int_t i=0; i<material_index; i++) {
      hamcTgtSlab *slab = components[i];
      rlen += slab->GetFracRadLen();
  }

  hamcTgtSlab *slab = components[material_index];
  zloff = slab->GetZloc() - 0.5*slab->GetLen();
  zdist = zscatt - zloff;
  rlen += slab->GetFracRadLen() * zdist / slab->GetLen();
  if (rlen < 0) cout << "hamcRadIn::ERROR: negative rlen ?"<<endl;

//   cout << "zlocation\t" << slab->GetZloc() << endl;
//   cout << "zscatt\t" << zscatt << endl;

  return rlen;

}   


Float_t hamcTarget::GetRadOut() {
// material_index = index of material where scattering occured
// Returns the fractional radiation length going out from scattering point.

  if (material_index == -1) return 0;
  Float_t rlen = 0;
  Float_t zloff = 0;
  Float_t zdist = 0;
  
  for (Int_t i=((Int_t)components.size()-1); i>material_index; i--) {
      hamcTgtSlab *slab = components[i];
      rlen += slab->GetFracRadLen();
  }

  hamcTgtSlab *slab = components[material_index];
  zloff = slab->GetZloc() + 0.5*slab->GetLen(); 
  zdist = zloff - zscatt;
  rlen += slab->GetFracRadLen() * zdist / slab->GetLen();
  if (rlen < 0) cout << "hamcRadOut::ERROR: negative rlen ?"<<endl;
  return rlen;

}   

Float_t hamcTarget::GetLenInMtl() {
// For the material that got struck, get the length into
// the material where it was hit
// "into" means along beam from front face to scatter point

  if (material_index == -1) return 0;
  hamcTgtSlab *slab = components[material_index];
  Float_t zloff = slab->GetZloc() - 0.5*slab->GetLen();
  Float_t zdist = zscatt - zloff;
  return zdist;

}


Float_t hamcTarget::GetLenOutMtl() {
// For the material that got struck, get the length out
// of the material where it was hit
// "out" means along beam from scatter point to exit of material

  if (material_index == -1) return 0;
  hamcTgtSlab *slab = components[material_index];
  Float_t zloff = slab->GetZloc() + 0.5*slab->GetLen(); 
  Float_t zdist = zloff - zscatt;
  return zdist;

}


Int_t hamcTarget::FindMtlIndex(Float_t zloc) {
// Find the index of the material for this Z location.
// The index tracks with (A,Z) and is use to lookup cross-section.
// The mass and atomic number are also looked up, as these
// are needed for kinematics and cross section, etc.

    Float_t zdiff;

    material_index = -1;
    mass = 0; 
    ascatt = -1;

    for (Int_t i=0; i<(Int_t)components.size(); i++) {

       hamcTgtSlab *slab = components[i];
       zdiff = zloc - slab->GetZloc();
       if (zdiff < 0) zdiff = -1.*zdiff;
// within 1/2 of length, the 0.502 is to account for roundoff
       if (zdiff < 0.502*slab->GetLen()) {
         material_index = i;
         ascatt = slab->GetA();
         mass = slab->GetMass();
         return OK;
       }

    }


// This would be a pretty serious error.  Should fix it if it happens
    cout << "hamcTarget::FindMtlIndex:WARNING: No component found at this Z " << zloc <<endl;
 
    return ERROR;

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

Float_t hamcTarget::GetMtlLen(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetLen();

}

Float_t hamcTarget::GetMtlEffLen(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetEffLen();

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

string hamcTarget::GetMtlName(Int_t iloc) {
 
  if (CheckIndex(iloc) == ERROR) return 0;

  return components[iloc]->GetName();

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





