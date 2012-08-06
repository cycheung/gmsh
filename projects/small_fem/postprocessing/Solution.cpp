#include "Solution.h"


using namespace std;

Solution::Solution(const System& system){
  // Get FunctionSpace
  const FunctionSpace& fs = system.getFunctionSpace();

  // Save some data
  this->sol     = &(system.getSol());
  this->dofM    = &(system.getDofManager());
  this->god     = &(dofM->getAllGroups());
  this->element = &(fs.getSupport().getAll());
  nGod          = god->size();
  nVertex       = fs.getSupport().getNVertex();

  // Interpolate
  switch(fs.getType()){
  case 0: 
  case 3: 
    nodalVectorValue = NULL;
    isScalar         = true;
    fsScalar         = 
      static_cast<const FunctionSpaceScalar*>(&fs);

    interpolateScalar();
    break;
    
  case 1: 
  case 2: 
    nodalScalarValue = NULL;
    isScalar         = false;

    interpolateVector();
    break;
    
  default: 
    // Impossible to be here
    break; 
  }
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
  // Alloc memory & set to zero
  nodalScalarValue = new vector<double>(nVertex);
  
  for(int i = 0; i < nVertex; i++)
    (*nodalScalarValue)[i] = 0;

  for(int i = 0; i < nGod; i++){
    // Get Dof
    const vector<Dof*>&  dof = (*god)[i]->getAll(); 
    const int           nDof = dof.size();

    // Get Coef
    vector<double> coef(nDof);
    for(int j = 0; j < nDof; j++)
      coef[j] = (*sol)(dofM->getGlobalId(*dof[j]));

    // Interpolate
    fsScalar->interpolateAtNodes((*god)[i]->getGeoElement(),
				 coef, 
				 *nodalScalarValue);
  }
}

void Solution::interpolateVector(void){
  // Alloc memory
  nodalVectorValue = new vector<fullVector<double> >(nVertex);
}
