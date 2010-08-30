//
// C++ Interface: partDomain
//
// Description: Class to replace an ElasticField (or a DGelasticField)
//              this class contains a bilinear and a linear terms for the iterative scheme
//              the data is also reorganized compared to elasticField
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PARTDOMAIN_H_
#define PARTDOMAIN_H_
#include "materialLaw.h"
#include "SVector3.h"
#include "DgC0PlateFunctionSpace.h"
#include "MInterfaceElement.h"
#include"Bindings.h"
#include "groupOfElements.h"
template<class T1, class T2> class DgC0BilinearTerm;
template<class T1> class DgC0LinearTerm;
class displacementField;
template<class T1,class T2> class IPField;
struct SolElementType{
  enum eltype{ShellPlaneStress, ShellPlaneStressWTI, ShellPlaneStressWF};
};

// class for a domain (pure virtual class)
class partDomain{
 protected :
  int _tag; // tag for the dofManager
  int _phys; // physical of interface group I don't know how to retrieve it from *g
  int _matnum; // number of materialLaw
  materialLaw *mat;
  DgC0BilinearTerm<SVector3,SVector3>* btermBulk;
  DgC0LinearTerm<SVector3>* ltermBulk;
  GaussQuadrature* integBulk;
  SolElementType::eltype _elemType;
  // bool interface term (allow to know if interface Terms have to be taken into account)
  bool interfaceterm;

 public :
  // Todo protect this variable
  groupOfElements *g; // support for this field

  // Constructors
  partDomain() : g(0), _tag(0), interfaceterm(false){}
  partDomain(const partDomain &source);
  partDomain& operator=(const partDomain &source){
    _tag = source._tag;
    _phys = source._phys;
    _matnum = source._matnum;
    btermBulk = source.btermBulk;
    ltermBulk = source.ltermBulk;
    integBulk = source.integBulk;
    _elemType = source._elemType;
    interfaceterm = source.interfaceterm;
    g = source.g;
    return *this;
  }
  //~partDomain(){delete btermBulk; delete ltermBulk; delete integBulk;}
  DgC0BilinearTerm<SVector3,SVector3>* getBilinearBulkTerm() const;
  DgC0LinearTerm<SVector3>* getLinearBulkTerm() const;
  GaussQuadrature* getBulkGaussIntegrationRule() const;
  SolElementType::eltype getSolElemType() const;
  int getTag()const;
  int getLawNum();
  int getPhysical() const;
  materialLaw* getMaterialLaw();
  bool IsInterfaceTerms() const;
  static void registerBindings(binding *b);
  virtual void setGaussIntegrationRule()=0;
  virtual void inverseLinearTermSign()=0;
// No virtual pure but lua problems ?? TODO fix it
  virtual void setTag(const int t)=0;
  virtual void setMaterialLaw(materialLaw* mlaw)=0;
  virtual void setLawNum(const int num)=0;
  virtual void setPhysical(const int ph)=0;
//
};
// class for Dg part domain (pure virtual)
class dgPartDomain : public partDomain{
 protected :
  DgC0BilinearTerm<SVector3,SVector3>* btermBound;
  DgC0LinearTerm<SVector3>* ltermBound;
  DgC0BilinearTerm<SVector3,SVector3>* btermVirtBound;
  DgC0LinearTerm<SVector3>* ltermVirtBound;
  GaussQuadrature* integBound;

 public :
  // TODO protect these variables
  std::vector<MInterfaceElement*> gi; // support for the interfaceElement TODO cast to a groupOfElements
  std::vector<MInterfaceElement*> gib; // support for the interfaceElement TODO cast to a groupOfElements

  dgPartDomain() : partDomain(){interfaceterm=true;}
  dgPartDomain(const dgPartDomain &source);
  dgPartDomain& operator=(dgPartDomain &source){
    this->partDomain::operator=(source);
    btermBound = source.btermBound;
    ltermBound = source.ltermBound;
    btermVirtBound = source.btermVirtBound;
    ltermVirtBound = source.ltermVirtBound;
    integBound = source.integBound;
    gi = source.gi;
    gib = source.gib;
    return *this;
  }
  //~dgPartDomain(){delete btermBound; delete ltermBound; delete btermVirtBound; delete ltermVirtBound; delete integBound;}
  DgC0BilinearTerm<SVector3,SVector3>* getBilinearInterfaceTerm();
  DgC0LinearTerm<SVector3>* getLinearInterfaceTerm()const;
  DgC0BilinearTerm<SVector3,SVector3>* getBilinearVirtualInterfaceTerm();
  DgC0LinearTerm<SVector3>* getLinearVirtualInterfaceTerm();
  GaussQuadrature* getInterfaceGaussIntegrationRule() const;
  virtual void setGaussIntegrationRule()=0;
  static void registerBindings(binding *b);
};

//#include "partDomain.h"
// Elasticfield for DG element (has a group of elements for interfaceelement)
class dgLinearShellDomain : public dgPartDomain{
 protected :
  double _h;
  int _msimp; // (odd) number of Simpson's point for integration on thickness
  bool FullDg; // To know which formulation Cg/Dg or FullDg is used for this field
 public :
// These functions have to be in part domain but binding problems with lua TODO fix it ??
  void setTag(const int t);
  void setMaterialLaw(materialLaw* mlaw);
  void setLawNum(const int num);
  void setPhysical(const int ph);
//
  void setFormulation(int b);
  bool getFormulation();
  int getmsimp()const;
  dgLinearShellDomain();
  dgLinearShellDomain(const dgLinearShellDomain &source);
  // problem with definition of operator=
  //~dgLinearShellDomain(){delete mat; delete g;} //If I delete mat here plante why ?? delete g ??
  double getThickness() const;
  virtual void setGaussIntegrationRule();
  void InitializeTerms(DgC0FunctionSpace<SVector3>& space1_,displacementField *uf,IPField<partDomain*,DgC0FunctionSpace<SVector3> >*ip,
                       double beta1, double beta2, double beta3);
  // lua interaction
  // function used by lua to set the properties
  void setThickness(const double h);
  void setNumberOfSimpsonPoints(int num);
  static void registerBindings(binding *b);
  void inverseLinearTermSign();
};


#endif // partdomain.h
