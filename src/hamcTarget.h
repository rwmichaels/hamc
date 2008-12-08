#ifndef ROOT_hamcTarget
#define ROOT_hamcTarget

//  hamcTarget   -- target in hamc
//  R. Michaels  June 2008

#include "Rtypes.h"
#include <vector>
#include <string>
#include <map>

class hamcExpt;

class hamcTgtSlab {
//  A monolithic slab of material belonging to the target
//  Private utility class used by hamcTarget
public:
// c'tor arguments are
// name,  index (tracks uniquely with [A,Z], needed by cross section)
// physical length 'len',  Normally in cm.
// effective length 'elen' includes estimated rad tail effect,
// radiation length 'rlen', atomic number 'a', charge 'z',
// mass 'm' (GeV), density 'd' (g/cm^3), Zlocation (cm)
  hamcTgtSlab(std::string mtl, Int_t id, Float_t len, Float_t elen, Float_t rlen, Float_t a, Float_t z, Float_t m, Float_t d): material(mtl),index(id),mtl_len(len), eff_len(elen),rad_len(rlen),anum(a),znum(z),mass(m),density(d),zlocation(0) { fracrl=0; if(rad_len!=0) fracrl=len/rad_len; };
  void  PutZloc(Float_t zz) { zlocation=zz; };
  Float_t GetZloc() const { return zlocation; };
  std::string GetName() const { return material; };
  Int_t   GetIndex() const { return index; };   // index, tracks with A,Z 
  Float_t GetLen() const { return mtl_len; };
  Float_t GetEffLen() const { return eff_len; };
  Float_t GetFracRadLen() const { return fracrl; };
  Float_t GetRadLen() const { return rad_len; };
  Float_t GetA() const { return anum; };
  Float_t GetZ() const { return znum; };
  Float_t GetMass() const { return mass; };
  Float_t GetDensity() const { return density; };
private:
  std::string material;
  Int_t index;
  Float_t mtl_len, eff_len, fracrl, rad_len, anum, znum, mass, density;
  Float_t zlocation;
};


class hamcTarget {

// A target.  Can be a composite of different materials.

  public:

     hamcTarget(std::string tgt_name);
     virtual ~hamcTarget()=0;       // abstract class

     virtual Int_t Init(hamcExpt*)=0;
// Generate a scattering location (do this each event)
     virtual Int_t Zscatt();      
// Z location of Scattering.
     Float_t GetZScatt() const { return zscatt;};   
// Material of scattering location (tracks with [A,Z] -- need for crsec)  
     Int_t   GetMtlIndex() const { return material_index; };
// Mass of scattering location (tracks with [A,Z] -- need for crsec)  
     Float_t GetMass() const { return mass; };

     virtual void Print();
     std::string GetName() const { return name; };
     Float_t GetLength() const { return overall_length; };
     Float_t GetEffLength() const { return effective_length; };
     Float_t GetA() const { return weighted_anum; };
     Float_t GetZ() const { return weighted_znum; };
     Float_t GetRadLength() const { return radiation_length; };
     Float_t GetDensity() const { return weighted_density; };
     Int_t   GetNumMtl() const { return components.size(); }; 
     Float_t GetMtlZloc(Int_t iloc);  // returns Z location of material #iloc
     Float_t GetMtlRadLen(Int_t iloc);
     Float_t GetMtlDensity(Int_t iloc);
     Float_t GetMtlA(Int_t iloc);
     Float_t GetMtlZ(Int_t iloc);
// atomic number of nucleus involved in scattering
     Float_t GetAscatt() { return ascatt; };
     Float_t GetRadIn();  // radiation length before scatt point
     Float_t GetRadOut(); // radiation length after scatt point.

  protected:

     std::string name;
     Bool_t did_init;
     Float_t overall_length, effective_length, radiation_length;
     Float_t width, weighted_anum, weighted_znum;
     Float_t weighted_density;
     std::vector<hamcTgtSlab *> components;
     Float_t zscatt, mass;
     Float_t ascatt;   
     Int_t material_index;

     Int_t Setup();
     Int_t CheckIndex(Int_t iloc);
     Int_t FindMtlIndex(Float_t z);  

  private: 

     hamcTarget& operator=(const hamcTarget& tgt);
     hamcTarget(const hamcTarget& tgt);


#ifndef NODICT
ClassDef (hamcTarget, 0)   // Target
#endif


};



#endif

   
