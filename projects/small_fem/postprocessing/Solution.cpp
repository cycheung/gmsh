#include "FunctionSpace.h"
#include "Dof.h"
#include "Solution.h"

using namespace std;

Solution::Solution(const System& system){
  // Get FunctionSpace
  const FunctionSpace& fs = system.getFunctionSpace();

  // Save some data
  this->sol     = &(system.getSol());
  this->dofM    = &(system.getDofManager());
  this->element = &(fs.getSupport().getAll());

  isScalar = true;
  interpolateScalar();

  for(int i = 0; i < nDof; i++)
    ;//printf("%lf\n", (*nodalScalarValue)[i]);
}

Solution::~Solution(void){
  if(isScalar)
    delete nodalScalarValue;

  else
    delete nodalVectorValue;
}

void Solution::write(const std::string name, Writer& writer) const{
  // Set Writer
  writer.setDomain(*element);

  if(isScalar)
    writer.setValues(*nodalScalarValue);

  else
    writer.setValues(*nodalVectorValue);

  // Write
  writer.write(name);
}

void Solution::interpolateScalar(void){
  const vector<const Dof*> dof = dofM->getAllDofs();
  nDof                         = dof.size();

  nodalScalarValue = new vector<double>(nDof);

  for(int i = 0; i < nDof; i++){
    (*nodalScalarValue)[dof[i]->getEntity() - 1] = 
      (*sol)(dofM->getGlobalId(*dof[i]));
  }
}
