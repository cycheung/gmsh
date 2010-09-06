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
#include "partDomain.h"
class GModelWithInterface;
class PView;
class groupOfElements;
class binding;
class displacementField;
class IPField;
template<class T1,class T2> class DgC0BilinearTerm;
template<class T1> class DgC0LinearTerm;

struct BoundaryCondition
{
  enum location{UNDEF,ON_VERTEX,ON_EDGE,ON_FACE,ON_VOLUME,PRESSURE};
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
  GaussQuadrature* integ;
  neumannBC () : BoundaryCondition(),_f(SVector3(0,0,0)){}
  void setGaussIntegrationRule(){
    groupOfElements::elementContainer::iterator it  = g->begin();
    MElement *ele = *it;
    if(ele->getPolynomialOrder()<3)
      integ = new GaussQuadrature(GaussQuadrature::Val);
    else
      integ = new GaussQuadrature(GaussQuadrature::ValVal);
  }
  //~neumannBC(){delete integ;} // segmentation fault si delete ?? car push_back in allNeumann puis "neu" est effacé lors de l'initialisation
};

class DgC0PlateSolver
{
 protected:
  GModelWithInterface *pModel;
  int _dim, _tag;
  dofManager<double> *pAssembler;
  FunctionSpace<SVector3> *LagSpace;
  dgGroupCollection _groups;
  std::vector<partDomain*> domainVector;
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
  bool hasInterface;
  bool fullDg; // Used to know which functionSpace has to be created. Remove this when FunctionSpace will be declared in partDomain

  // Function to initiate the solver before solve
  // (create interfaceElement and link the materialLaw and dgLinearShellDomain)
  void init();

  // Function used by non linear solver
  void NewtonRaphson(linearSystem<double> *lsys, dofManager<double> *pAssembler, displacementField *ufield,
                     IPField *ipf,
                     std::vector<MInterfaceElement*> &vinter);

  double computeNorm0(linearSystem<double> *lsys, dofManager<double> *pAssembler, displacementField *ufield,
                      IPField *ipf,
                      std::vector<MInterfaceElement*> &vinter);

  double computeRightHandSide(linearSystem<double> *lsys, dofManager<double> *pAssembler, displacementField *ufield,
                                    IPField *ipf,
                                    std::vector<MInterfaceElement*> &vinter);

  void computeStiffMatrix(linearSystem<double> *lsys, dofManager<double> *pAssembler,
                                    displacementField *ufield, IPField *ipf);

 public : //protected with lua ??
  enum solver{ Gmm=0,Taucs=1,Petsc=2};
  enum scheme{StaticLinear=0, StaticNonLinear=1};
  solver whatSolver; //  Solver used to solve
  scheme whatScheme; // scheme used to solve equation
 public:
  DgC0PlateSolver(int tag) : _tag(tag), fullDg(false), LagSpace(0),pAssembler(0), numstep(1), endtime(1.), _tol(1.e-6), nsba(1),
                             whatSolver(DgC0PlateSolver::Gmm), whatScheme(DgC0PlateSolver::StaticLinear),
                             _beta1(10.), _beta2(10.), _beta3(10.), hasInterface(false){} // default Gmm and static linear scheme

  virtual ~DgC0PlateSolver()
  {
    if (LagSpace) delete LagSpace;
    if (pAssembler) delete pAssembler;
    for(std::map<int,materialLaw*>::iterator it = maplaw.begin(); it!=maplaw.end();++it){ delete it->second;}
    for(std::vector<partDomain*>::iterator it = domainVector.begin(); it!=domainVector.end(); ++it){delete *it;}
  }
  void readInputFile(const std::string meshFileName);
  void setMesh(const std::string meshFileName);
  void setSolver(const int s){whatSolver= (solver)s;}
  void setScheme(const int s){whatScheme=(scheme)s;}
  void setSNLData(const int ns, const double et, const double reltol);
  void setStepBetweenArchiving(const int n){nsba = n;}
  void setFormulation(const int f){ f==0 ? fullDg = false : fullDg = true;}
  bool getFormulation() const{return fullDg;}
  int getStepBetweenArchiving() const{return nsba;}
  int getNumStep() const{return numstep;}
  double getEndTime() const{return endtime;}
  scheme getScheme() const{return whatScheme;}
  virtual void solve();
  virtual void solveSNL();
  virtual void setStabilityParameters(const double b1, const double b2=0., const double b3=0.){_beta1=b1;_beta2=b2;_beta3=b3;}
  virtual void createInterfaceElement();
  virtual materialLaw* getMaterialLaw(const int num);
  virtual bool IsInterface() const{return hasInterface;}

  // create interfaceelement with dgGoupOfElement from dg project doesn't work (segmentation fault)
  virtual void createInterfaceElement_2();

  virtual std::vector<Dof> getDofArchForce();
  // functions for lua interaction
  virtual void solveLUA();
  virtual void addDgLinearElasticShellDomain(dgLinearShellDomain* ED, const int f, const int d);
  virtual void addMaterialLaw(materialLaw *mlaw);
  virtual void addLinearElasticLawPlaneStress(linearElasticLawPlaneStress *mlaw);
  virtual void addLinearElasticLawPlaneStressWithFracture(linearElasticLawPlaneStressWithFracture *mlaw);
  virtual void addTheta(const int numphys);
  virtual void addDisp(std::string onwhat, const int numphys, const int comp, const double value);
  virtual void addForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval);
  virtual void addIndepDisp(std::string onwhat, const int numphys, const int comp, const double value);
  virtual void addIndepForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval);
  virtual void addArchivingForceForPhysicalGroup(const int numphys, const int comp);
  virtual void addArchivingNodeDisplacement(const int num, const int comp);
  virtual void addPhysInitBroken(const int phys);
  virtual void addPressureOnPhysicalGroup(const int numphys, const double press);
  static void registerBindings(binding *b);
};


#endif
