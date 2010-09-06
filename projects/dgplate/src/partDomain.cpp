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
#include "DgC0PlateSolver.h"
#include "DgC0PlateTerms.h"
partDomain::partDomain(const partDomain &source){
  _tag = source._tag;
  _phys = source._phys;
  _matnum = source._matnum;
  btermBulk = source.btermBulk;
  ltermBulk = source.ltermBulk;
  integBulk = source.integBulk;
  _elemType = source._elemType;
  interfaceterm = source.interfaceterm;
  g = source.g;
}

DgC0BilinearTerm<SVector3,SVector3>* partDomain::getBilinearBulkTerm() const{return btermBulk;}
DgC0LinearTerm<SVector3>* partDomain::getLinearBulkTerm() const{return ltermBulk;}
GaussQuadrature* partDomain::getBulkGaussIntegrationRule() const{return integBulk;}
SolElementType::eltype partDomain::getSolElemType() const {return _elemType;}
bool partDomain::IsInterfaceTerms() const{return interfaceterm;}

int partDomain::getTag()const{return _tag;}
int partDomain::getLawNum(){return _matnum;}
int partDomain::getPhysical() const{return _phys;}
materialLaw* partDomain::getMaterialLaw() {return mat;}

dgPartDomain::dgPartDomain(const dgPartDomain &source) : partDomain(source){
  btermBound = source.btermBound;
  ltermBound = source.ltermBound;
  btermVirtBound = source.btermVirtBound;
  ltermVirtBound = source.ltermVirtBound;
  integBound = source.integBound;
  gi = source.gi;
  gib = source.gib;
}
DgC0BilinearTerm<SVector3,SVector3>* dgPartDomain::getBilinearInterfaceTerm(){return btermBound;}
DgC0LinearTerm<SVector3>* dgPartDomain::getLinearInterfaceTerm()const{return ltermBound;}
DgC0BilinearTerm<SVector3,SVector3>* dgPartDomain::getBilinearVirtualInterfaceTerm(){return btermVirtBound;}
DgC0LinearTerm<SVector3>* dgPartDomain::getLinearVirtualInterfaceTerm(){return ltermVirtBound;}
GaussQuadrature* dgPartDomain::getInterfaceGaussIntegrationRule() const {return integBound;}


bool dgLinearShellDomain::getFormulation(){return FullDg;}
int dgLinearShellDomain::getmsimp()const {return _msimp;}
dgLinearShellDomain::dgLinearShellDomain () : dgPartDomain(){}
dgLinearShellDomain::dgLinearShellDomain(const dgLinearShellDomain &source) : dgPartDomain(source){
  _h = source._h;
  _msimp = source._msimp;
  FullDg = source.FullDg;
}

void dgLinearShellDomain::setFormulation(int b){b==1 ? FullDg=true : FullDg=false;}
double dgLinearShellDomain::getThickness() const{return _h;}
void dgLinearShellDomain::setThickness(const double h){_h=h;}
void dgLinearShellDomain::setNumberOfSimpsonPoints(int num){num%2==0 ? _msimp=++num : _msimp=num;}
void dgLinearShellDomain::setMaterialLaw(materialLaw *mlaw){
  mat = mlaw;
  switch(mat->getType()){
    case materialLaw::linearElasticPlaneStress :
      if(_msimp==1) _elemType = SolElementType::ShellPlaneStress; // change this
      else _elemType = SolElementType::ShellPlaneStressWTI;
      break;
    case materialLaw::linearElasticPlaneStressWithFracture :
      _elemType = SolElementType::ShellPlaneStressWF;
      if(_msimp==1){
        _msimp=3;
        Msg::Warning("Impossible to compute a fracture problem with only 1 Simpson's point. The number of Simpson's point is put to 3.\n");
      }
      if(!FullDg){
        FullDg = true;
        Msg::Warning("Impossible to compute a fracture problem with cG/dG formulation. Switch to full dG formulation\n");
      }
      break;
  }
}

void dgLinearShellDomain::InitializeTerms(FunctionSpace<SVector3>& space1_,displacementField *uf,
                                         IPField*ip,
                                         double beta1, double beta2, double beta3){
  DgC0FunctionSpace<SVector3>* dgspace = dynamic_cast<DgC0FunctionSpace<SVector3>*>(&space1_);
  btermBulk = new IsotropicElasticStiffBulkTermC0Plate(*dgspace,mat,_h,FullDg,uf,ip,_elemType);
  btermBound = new IsotropicElasticStiffInterfaceTermC0Plate(*dgspace,mat,beta1,beta2,beta3, _h,uf,ip,_elemType,FullDg,1.e-8);
  // remove the true ??
  btermVirtBound = new IsotropicElasticStiffVirtualInterfaceTermC0Plate(*dgspace,mat,beta1,beta2,beta3,_h,uf,ip,_elemType, FullDg);
  ltermBulk = new IsotropicElasticForceBulkTermC0Plate(*dgspace,mat,_h,FullDg,uf,ip,_elemType);
  ltermBound = new IsotropicElasticForceInterfaceTermC0Plate(*dgspace, mat,beta1,beta2,beta3,_h,FullDg,uf,ip,_elemType);
  ltermVirtBound = new IsotropicElasticForceVirtualInterfaceTermC0Plate(*dgspace,mat,beta1,beta2,beta3,_h, FullDg,uf,ip,_elemType);
}

void dgLinearShellDomain::setLawNum(const int num){_matnum = num;}

void partDomain::registerBindings(binding *b){
  classBinding *cb = b->addClass<partDomain>("partDomain"); //TODO create MaterialField and derive dgLinearShellDomain
  cb->setDescription("part Domain");
}
void dgPartDomain::registerBindings(binding *b){
  classBinding *cb = b->addClass<dgPartDomain>("dgPartDomain"); //TODO create MaterialField and derive dgLinearShellDomain
  cb->setDescription("discontinuous galerkin part domain");
}

void dgLinearShellDomain::setGaussIntegrationRule(){
  groupOfElements::elementContainer::iterator it = g->begin();
  MElement *e = *it;
  if(e->getPolynomialOrder()<3){
    integBulk = new GaussQuadrature(GaussQuadrature::Grad);
    integBound = new GaussQuadrature(GaussQuadrature::Val);
  }
  else{
    integBulk = new GaussQuadrature(GaussQuadrature::GradGrad);
    integBound = new GaussQuadrature(GaussQuadrature::ValVal);
  }
}

void dgLinearShellDomain::inverseLinearTermSign(){
  ltermBulk->invSign();
  ltermBound->invSign();
  ltermVirtBound->invSign();
}

void dgLinearShellDomain::setTag(const int t){_tag=t;}
void dgLinearShellDomain::setMaterialLaw(materialLaw* mlaw);
void dgLinearShellDomain::setLawNum(const int num);
void dgLinearShellDomain::setPhysical(const int ph){_phys = ph;}

void dgLinearShellDomain::registerBindings(binding *b)
{
  classBinding *cb = b->addClass<dgLinearShellDomain>("dgLinearShellDomain"); //TODO create MaterialField and derive dgLinearShellDomain
  cb->setDescription("Material field: Give all parameters for material law and numerical integration");

  methodBinding *cm;
  // Constructor
  cm = cb->setConstructor<dgLinearShellDomain>();
  cm->setDescription("Initialisation of a linear elastic law. All the methods need to be used to set correctly this field");
  // method
  cm = cb->addMethod("thickness", &dgLinearShellDomain::setThickness);
  cm->setArgNames("h",NULL);
  cm->setDescription("Value of thickness (because beam, plate or shell formulation");
  cm = cb->addMethod("simpsonPoints", &dgLinearShellDomain::setNumberOfSimpsonPoints);
  cm->setArgNames("num",NULL);
  cm->setDescription("Set the number of Simpson's point used for this field. If = 1 no integration (linear elastic only). If the given is even 1 is add to this number");
  cm = cb->addMethod("tag", &dgLinearShellDomain::setTag);
  cm->setArgNames("t",NULL);
  cm->setDescription("Set the tag of this field");
  cm = cb->addMethod("formulation", &dgLinearShellDomain::setFormulation);
  cm->setArgNames("b",NULL);
  cm->setDescription("Set the material law used for this field. LinearElasticPlaneStress=0");
  cm = cb->addMethod("lawnumber", &dgLinearShellDomain::setLawNum);
  cm->setArgNames("num",NULL);
  cm->setDescription("Give the number of the law linked to this elasticField");
}
