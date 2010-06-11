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

// Elasticfield for DG element (has a group of elements for interfaceelement) USED with materialLaw change name and rewritten ??
struct DGelasticField{
  int _tag; // tag for the dofManager
  int _phys; // physical of interface group I don't know how to retrieve it from *g
  int _msimp; // (odd) number of Simpson's point for integration on thickness
  materialLaw::matname _matname; // material law used (elastic for now change class name) USED ??
  materialLaw *mat;
  groupOfElements *g; // support for this field
  std::vector<MInterfaceElement*> gi; // support for the interfaceElement TODO cast to a groupOfElements
  std::vector<MInterfaceElement*> gib; // support for the interfaceElement TODO cast to a groupOfElements
  SolElementType::eltype _elemType;
  double _E, _nu, _beta1, _beta2, _beta3, _h;  // specific elastic datas (should be somewhere else)
  double _Gc, _sigmac; // fracture parameter
  bool FullDg; // To know which formulation Cg/Dg or FullDg is used for this field
  void setFormulation(int b){b==1 ? FullDg=true : FullDg=false;}
  void setMaterialLaw(int m){
    _matname=materialLaw::matname(m);
    switch(m){
      case materialLaw::linearElasticPlaneStress :
        mat = new linearElasticLawPlaneStress(_E,_nu);
        if(_msimp==1) _elemType = SolElementType::PlatePlaneStress; // change this
        else _elemType = SolElementType::PlatePlaneStressWTI;
        break;
      case materialLaw::linearElasticPlaneStressWithFracture :
        mat = new linearElasticLawPlaneStressWithFracture(_E,_nu,_Gc,_sigmac);
        _elemType = SolElementType::PlatePlaneStressWF;
        if(_msimp==1){
          _msimp=3;
          Msg::Warning("Impossible to compute a fracture problem with only 1 Simpson's point. The number of Simpson's point is put to 3.\n");
        }
        if(!FullDg){
          FullDg = true;
          Msg::Warning("Impossible to compute a fracture problem with cG/dG formulation. Switch to full dG formulation\n");
        }
        break;
      default : printf("The chosen material law seems to be inexsistent\n");
    }
  }
  materialLaw* getMaterialLaw() {return mat;}
  bool getFormulation(){return FullDg;}
  double getSigmaC() const {return _sigmac;}
  SolElementType::eltype getSolElemType() const {return _elemType;}
  int getmsimp()const {return _msimp;}
  DGelasticField () : g(0), _tag(0){}
  //~DGelasticField(){delete mat; delete g;} //If I delete mat here plante why ?? delete g ??
  // lua interaction
  // function used by lua to set the properties
  void setYoungModulus(const double E){_E=E;}
  void setPoissonRation(const double nu){_nu=nu;}
  void setStabilityParameters(const double b1, const double b2=0., const double b3=0.){_beta1=b1;_beta2=b2;_beta3=b3;}
  void setThickness(const double h){_h=h;}
  void setNumberOfSimpsonPoints(int num){num%2==0 ? _msimp=++num : _msimp=num;}
  void setTag(const int t){_tag=t;}
  //void setElementGroup(const int phys, const int dim){_phys=phys; g = new groupOfElements(dim,phys);}
  void setFractureParameter(const double Gc_, const double sigmac_){_Gc = Gc_; _sigmac = sigmac_;}

  static void registerBindings(binding *b)
  {
  classBinding *cb = b->addClass<DGelasticField>("DGelasticField"); //TODO create MaterialField and derive DGelasticField
  cb->setDescription("Material field: Give all parameters for material law and numerical integration");

  methodBinding *cm;
  // Constructor
  cm = cb->setConstructor<DGelasticField>();
  cm->setDescription("Initialisation of a linear elastic law. All the methods need to be used to set correctly this field");
  // method
  cm = cb->addMethod("young", &DGelasticField::setYoungModulus);
  cm->setArgNames("E",NULL);
  cm->setDescription("Set the Young modulus value");
  cm = cb->addMethod("poisson", &DGelasticField::setPoissonRation);
  cm->setArgNames("nu",NULL);
  cm->setDescription("Set the Poisson ratio value");
  cm = cb->addMethod("stabilityParameters", &DGelasticField::setStabilityParameters);
  cm->setArgNames("b1","b2","b3",NULL);
  cm->setDescription("Beta value for stabilization. The first one only for Cg/Dg formulation");
  cm = cb->addMethod("thickness", &DGelasticField::setThickness);
  cm->setArgNames("h",NULL);
  cm->setDescription("Value of thickness (because beam, plate or shell formulation");
  cm = cb->addMethod("simpsonPoints", &DGelasticField::setNumberOfSimpsonPoints);
  cm->setArgNames("num",NULL);
  cm->setDescription("Set the number of Simpson's point used for this field. If = 1 no integration (linear elastic only). If the given is even 1 is add to this number");
  cm = cb->addMethod("tag", &DGelasticField::setTag);
  cm->setArgNames("t",NULL);
  cm->setDescription("Set the tag of this field");
  cm = cb->addMethod("formulation", &DGelasticField::setFormulation);
  cm->setArgNames("b",NULL);
  cm->setDescription("Int to set the formulation used for this field. Cg/Dg=0 full Dg=1");
  cm = cb->addMethod("law", &DGelasticField::setMaterialLaw);
  cm->setArgNames("m",NULL);
  cm->setDescription("Set the material law used for this field. LinearElasticPlaneStress=0");
  cm = cb->addMethod("fractureParameter", &DGelasticField::setFractureParameter);
  cm->setArgNames("Gc_","sigmac_",NULL);
  cm->setDescription("Set the fracture parameter. First argument Gc second argument sigma_c");
  }
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
  // map to archive force
  std::map<int,std::vector<Dof> > aef;
  std::map<int,double> aefvalue;

  // std vector to archive a node displacement
  std::vector<Dof> anoded;

  int numstep; // Number of step not used for StaticLinearScheme but it is not necesary to derive class because (small useless data)
  double endtime; // final time not used for StaticLinearScheme but it is not necesary to derive class because (small useless data)
  double _tol; // relative tolerance for iteration not used for StaticLinearScheme but it is not necesary to derive class because (small useless data)
  int nsba; // number of step between 2 archiving (=1 by default)

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
  DgC0PlateSolver(int tag) : _tag(tag),LagSpace(0),pAssembler(0), numstep(1), endtime(1.), _tol(1.e-6), nsba(1)
  {
    whatSolver=DgC0PlateSolver::Gmm; whatScheme=DgC0PlateSolver::StaticLinear; // default Gmm and static linear scheme
  }
  virtual ~DgC0PlateSolver()
  {
    if (LagSpace) delete LagSpace;
    if (pAssembler) delete pAssembler;
  }
  void readInputFile(const std::string meshFileName);
  void setMesh(const std::string meshFileName);
  void setSolver(const int s){whatSolver= (solver)s;}
  void setScheme(const int s){whatScheme=(scheme)s;}
  void setSNLData(const int ns, const double et, const double reltol){
    if(whatScheme==StaticNonLinear){
      numstep=ns;endtime=et;_tol=reltol;
    }
    else{
      Msg::Error("Impossible to set data for Static Non Linear scheme because another is chosen to solve the problem");
    }
  }
  void setStepBetweenArchiving(const int n){nsba = n;}
  int getStepBetweenArchiving() const{return nsba;}
  int getNumStep() const{return numstep;}
  double getEndTime() const{return endtime;}
  scheme getScheme() const{return whatScheme;}
  virtual void solve();
  virtual void solveSNL();
  virtual void createInterfaceElement();
  // create interfaceelement with dgGoupOfElement from dg projet doesn't work (segmentation fault)
  virtual void createInterfaceElement_2();

  virtual std::vector<Dof> getDofArchForce();
  // functions for lua interaction
  virtual void solveLUA(){
    switch(whatScheme){
      case StaticLinear :
        this->solve();
        break;
      case StaticNonLinear :
        this->solveSNL();
        break;
    }
  }
  virtual void addElasticDomain(DGelasticField* ED, const int f, const int d);
  virtual void addTheta(const int numphys);
  virtual void addDisp(std::string onwhat, const int numphys, const int comp, const double value){
    dirichletBC diri;
    const std::string node("Node");
    const std::string edge("Edge");
    const std::string face("Face");
    if(onwhat==node){
      diri.g = new groupOfElements (0, numphys);
      diri.onWhat=BoundaryCondition::ON_VERTEX;
    }
    else if(onwhat==edge){
      diri.g = new groupOfElements (1, numphys);
      diri.onWhat=BoundaryCondition::ON_EDGE;
    }
    else if(onwhat==face){
      diri.g = new groupOfElements (2, numphys);
      diri.onWhat=BoundaryCondition::ON_FACE;
    }
    else Msg::Error("Impossible to prescribe a displacement on a %s\n",onwhat.c_str());
    diri._f= simpleFunctionTime<double>(value);
    diri._comp=comp;
    diri._tag=numphys;
    allDirichlet.push_back(diri);
  }

  virtual void addForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval){
    neumannBC neu;
    const std::string node("Node");
    const std::string edge("Edge");
    const std::string face("Face");
    const std::string volume("Volume");
    if(onwhat==node){
      neu.g = new groupOfElements (0, numphys);
      neu.onWhat=BoundaryCondition::ON_VERTEX;
    }
    else if(onwhat==edge){
      neu.g = new groupOfElements (1, numphys);
      neu.onWhat=BoundaryCondition::ON_EDGE;
    }
    else if(onwhat==face){
      neu.g = new groupOfElements (2, numphys);
      neu.onWhat=BoundaryCondition::ON_FACE;
    }
    else if(onwhat==volume){
      neu.g = new groupOfElements (3, numphys);
      neu.onWhat=BoundaryCondition::ON_VOLUME;
    }
    else  Msg::Error("Impossible to prescribe a force on a %s\n",onwhat.c_str());
    neu._f= simpleFunctionTime<SVector3>(SVector3(xval, yval, zval));
    neu._tag=numphys;
    allNeumann.push_back(neu);
  }
  virtual void addArchivingEdgeForce(const int numphys, const int comp);
  virtual void addArchivingNodeDisplacement(const int num, const int comp);
  static void registerBindings(binding *b);
};


#endif
