// Gmsh - Copyright (C) 1997-2009 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _DGC0PLATE_SOLVER_H_
#define _DGC0PLATE_SOLVER_H_

#include <map>
#include <string>
#include "SVector3.h"
#include "dofManager.h"
#include "simpleFunction.h"
#include "DgC0PlateFunctionSpace.h"
#include "MInterfaceElement.h"
#include "../dg/dgGroupOfElements.h"
#include "materialLaw.h"

class GModelWithInterface;
class PView;
class groupOfElements;
class binding;
struct elasticField {
  int _tag; // tag for the dofManager
  groupOfElements *g; // support for this field
  double _E, _nu; // specific elastic datas (should be somewhere else)
  elasticField () : g(0),_tag(0){}
};

struct SolElementType{
  enum eltype{PlatePlaneStress};
};

// Elasticfield for DG element (has a group of elements for interfaceelement) USED with materialLaw change name and rewritten ??
struct DGelasticField{
  int _tag; // tag for the dofManager
  int _phys; // physical of interface group I don't know how to retrieve it from *g
  materialLaw::matname _matname; // material law used (elastic for now change class name) USED ??
  materialLaw *mat;
  groupOfElements *g; // support for this field
  std::vector<MInterfaceElement*> gi; // support for the interfaceElement TODO cast to a groupOfElements
  std::vector<MInterfaceElement*> gib; // support for the interfaceElement TODO cast to a groupOfElements
  SolElementType::eltype _elemType;
  double _E, _nu, _beta1, _beta2, _beta3, _h;  // specific elastic datas (should be somewhere else)
  bool FullDg; // To know which formulation Cg/Dg or FullDg is used for this field
  void setFormulation(int b){b==1 ? FullDg=true : FullDg=false;}
  void setMaterialLaw(const materialLaw::matname m){
    _matname=m;
    switch(m){
      case materialLaw::linearElasticPlaneStress :
        mat = new linearElasticLawPlaneStress(_E,_nu);
        _elemType = SolElementType::PlatePlaneStress; // change this
        break;
      default : printf("The chosen material law seems to be inexsistent\n");
    }
  }
  materialLaw* getMaterialLaw() {return mat;}
  bool getFormulation(){return FullDg;}
  SolElementType::eltype getSolElemType(){return _elemType;}
  DGelasticField () : g(0), _tag(0){}
  //~DGelasticField(){delete mat;} //If I delete mat here plante why ??
};

struct BoundaryCondition
{
  enum location{UNDEF,ON_VERTEX,ON_EDGE,ON_FACE,ON_VOLUME};
  location onWhat; // on vertices or elements
  int _tag; // tag for the dofManager
  groupOfElements *g; // support for this BC
  BoundaryCondition() : g(0),_tag(0),onWhat(UNDEF) {};
};

struct dirichletBC : public BoundaryCondition
{
  int _comp; // component
  simpleFunction<double> _f;
  dirichletBC () :BoundaryCondition(),_comp(0),_f(0){}
};

struct neumannBC  : public BoundaryCondition
{
  simpleFunction<SVector3> _f;
  neumannBC () : BoundaryCondition(),_f(SVector3(0,0,0)){}
};

// an elastic solver ...
class DgC0PlateSolver
{
 protected:
  GModelWithInterface *pModel;
  int _dim, _tag;
  dofManager<double> *pAssembler;
  DgC0FunctionSpace<SVector3> *LagSpace;
  dgGroupCollection _groups;
  // young modulus and poisson coefficient per physical
  std::vector<DGelasticField> elasticFields;
  // neumann BC
  std::vector<neumannBC> allNeumann;
  // dirichlet BC
  std::vector<dirichletBC> allDirichlet;
  // typedef for IPvariable (should be somewhere else)

 public : //protected with lua ??
  enum solver{ Gmm=0,Taucs=1,Petsc=2};
  solver whatSolver; //  Solver used to solve

 public:
  DgC0PlateSolver(int tag) : _tag(tag),LagSpace(0),pAssembler(0) {whatSolver=DgC0PlateSolver::Gmm;} // default Gmm
  virtual ~DgC0PlateSolver()
  {
    if (LagSpace) delete LagSpace;
    if (pAssembler) delete pAssembler;
  }
  void readInputFile(const std::string &meshFileName);
  void setMesh(const std::string &meshFileName);
  void setSolver(const int s){whatSolver= (solver)s;}
  virtual void solve();
  virtual void createInterfaceElement();
  // create interfaceelement with dgGoupOfElement from dg projet doesn't work (segmentation fault)
  virtual void createInterfaceElement_2();
  virtual PView *buildDisplacementView(const std::string &postFileName);
  void buildVonMisesView(const std::string &postFileName);
  //virtual PView *buildElasticEnergyView(const std::string &postFileName);
  //virtual PView *buildVonMisesView(const std::string &postFileName);
  // std::pair<PView *, PView*> buildErrorEstimateView
  //   (const std::string &errorFileName, double, int);
  // std::pair<PView *, PView*> buildErrorEstimateView
  //   (const std::string &errorFileName, const elasticityData &ref, double, int);
  void registerBindings(binding *b);
};


#endif
