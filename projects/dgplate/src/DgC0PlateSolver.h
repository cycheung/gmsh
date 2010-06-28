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
#include"Bindings.h"
#include "SimpleFunctionTime.h"
#include "groupOfElements.h"
class GModelWithInterface;
class PView;
class groupOfElements;
class binding;
class displacementField;
template<class T1,class T2> class IPField;

struct elasticField {
  int _tag; // tag for the dofManager
  groupOfElements *g; // support for this field
  double _E, _nu; // specific elastic datas (should be somewhere else)
  elasticField () : g(0),_tag(0){}
};

struct SolElementType{
  enum eltype{PlatePlaneStress, PlatePlaneStressWTI, PlatePlaneStressWF};
};

//#include "partDomain.h"
// Elasticfield for DG element (has a group of elements for interfaceelement) USED with materialLaw change name and rewritten ??
struct DGelasticField{
  int _tag; // tag for the dofManager
  int _phys; // physical of interface group I don't know how to retrieve it from *g
  int _msimp; // (odd) number of Simpson's point for integration on thickness
  int _matnum; // number of materialLaw
  materialLaw *mat;
  groupOfElements *g; // support for this field
  std::vector<MInterfaceElement*> gi; // support for the interfaceElement TODO cast to a groupOfElements
  std::vector<MInterfaceElement*> gib; // support for the interfaceElement TODO cast to a groupOfElements
  SolElementType::eltype _elemType;
  double _h;  // specific elastic datas (should be somewhere else)
  bool FullDg; // To know which formulation Cg/Dg or FullDg is used for this field
  void setFormulation(int b);
//  void setMaterialLaw(int m);
  void setMaterialLaw(materialLaw* mlaw);
  int getLawNum(){return _matnum;}
  void setLawNum(const int num);
  materialLaw* getMaterialLaw() {return mat;}
  bool getFormulation(){return FullDg;}
  SolElementType::eltype getSolElemType() const {return _elemType;}
  int getmsimp()const {return _msimp;}
  DGelasticField () : g(0), _tag(0){}
  //~DGelasticField(){delete mat; delete g;} //If I delete mat here plante why ?? delete g ??

  // lua interaction
  // function used by lua to set the properties
  void setThickness(const double h){_h=h;}
  void setNumberOfSimpsonPoints(int num){num%2==0 ? _msimp=++num : _msimp=num;}
  void setTag(const int t){_tag=t;}
  static void registerBindings(binding *b);
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
  simpleFunctionTime<double> _f;
  dirichletBC () :BoundaryCondition(),_comp(0),_f(0){}
};

struct neumannBC  : public BoundaryCondition
{
  simpleFunctionTime<SVector3> _f;
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
  // vector with material law
  std::map<int,materialLaw*> maplaw;
  // physical entities that are initialy broken
  std::vector<int> initbrokeninter;


  // map to archive force
  std::map<int,std::vector<Dof> > aef;
  std::map<int,double> aefvalue;
  // std vector to archive a node displacement
  std::vector<Dof> anoded;

  // specific data
  double _beta1, _beta2, _beta3; // Stability parameters for cG/dG case only beta1 is used
  int numstep; // Number of step not used for StaticLinearScheme but it is not necesary to derive class because (small useless data)
  double endtime; // final time not used for StaticLinearScheme but it is not necesary to derive class because (small useless data)
  double _tol; // relative tolerance for iteration not used for StaticLinearScheme but it is not necesary to derive class because (small useless data)
  int nsba; // number of step between 2 archiving (=1 by default)

  // Function to initiate the solver before solve
  // (create interfaceElement and link the materialLaw and DGElasticField)
  void init();

  // Function used by non linear solver
  void NewtonRaphson(linearSystem<double> *lsys, dofManager<double> *pAssembler, displacementField *ufield,
                     IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf, GaussQuadrature &integbound, GaussQuadrature &integbulk,
                     std::vector<MInterfaceElement*> &vinter);

  double computeNorm0(linearSystem<double> *lsys, dofManager<double> *pAssembler, displacementField *ufield,
                      IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf,
                      GaussQuadrature &integbound, GaussQuadrature &integbulk, std::vector<MInterfaceElement*> &vinter);

  double computeRightHandSide(linearSystem<double> *lsys, dofManager<double> *pAssembler, displacementField *ufield,
                                    IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf, GaussQuadrature &integbound, GaussQuadrature &integbulk,
                                    std::vector<MInterfaceElement*> &vinter);

  void computeStiffMatrix(linearSystem<double> *lsys, dofManager<double> *pAssembler,
                                    displacementField *ufield, IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf,
                                    GaussQuadrature &Integ_Boundary, GaussQuadrature &Integ_Bulk,
                                    std::vector<MInterfaceElement*> &vinter);

 public : //protected with lua ??
  enum solver{ Gmm=0,Taucs=1,Petsc=2};
  enum scheme{StaticLinear=0, StaticNonLinear=1};
  solver whatSolver; //  Solver used to solve
  scheme whatScheme; // scheme used to solve equation
 public:
  DgC0PlateSolver(int tag) : _tag(tag),LagSpace(0),pAssembler(0), numstep(1), endtime(1.), _tol(1.e-6), nsba(1),
                             whatSolver(DgC0PlateSolver::Gmm), whatScheme(DgC0PlateSolver::StaticLinear),
                             _beta1(10.), _beta2(10.), _beta3(10.){} // default Gmm and static linear scheme

  virtual ~DgC0PlateSolver()
  {
    if (LagSpace) delete LagSpace;
    if (pAssembler) delete pAssembler;
    for(std::map<int,materialLaw*>::iterator it = maplaw.begin(); it!=maplaw.end();++it){ delete it->second;}
  }
  void readInputFile(const std::string meshFileName);
  void setMesh(const std::string meshFileName);
  void setSolver(const int s){whatSolver= (solver)s;}
  void setScheme(const int s){whatScheme=(scheme)s;}
  void setSNLData(const int ns, const double et, const double reltol);
  void setStepBetweenArchiving(const int n){nsba = n;}
  int getStepBetweenArchiving() const{return nsba;}
  int getNumStep() const{return numstep;}
  double getEndTime() const{return endtime;}
  scheme getScheme() const{return whatScheme;}
  virtual void solve();
  virtual void solveSNL();
  virtual void setStabilityParameters(const double b1, const double b2=0., const double b3=0.){_beta1=b1;_beta2=b2;_beta3=b3;}
  virtual void createInterfaceElement();
  virtual materialLaw* getMaterialLaw(const int num);
  // create interfaceelement with dgGoupOfElement from dg projet doesn't work (segmentation fault)
  virtual void createInterfaceElement_2();

  virtual std::vector<Dof> getDofArchForce();
  // functions for lua interaction
  virtual void solveLUA();
  virtual void addElasticDomain(DGelasticField* ED, const int f, const int d);
  virtual void addMaterialLaw(materialLaw *mlaw);
  virtual void addLinearElasticLawPlaneStress(linearElasticLawPlaneStress *mlaw);
  virtual void addLinearElasticLawPlaneStressWithFracture(linearElasticLawPlaneStressWithFracture *mlaw);
  virtual void addTheta(const int numphys);
  virtual void addDisp(std::string onwhat, const int numphys, const int comp, const double value);
  virtual void addForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval);
  virtual void addIndepDisp(std::string onwhat, const int numphys, const int comp, const double value);
  virtual void addIndepForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval);
  virtual void addArchivingEdgeForce(const int numphys, const int comp);
  virtual void addArchivingNodeDisplacement(const int num, const int comp);
  virtual void addPhysInitBroken(const int phys);
  static void registerBindings(binding *b);
};


#endif
