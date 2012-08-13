#include "Solution.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "MElement.h"
#include "MVertex.h"
#include "Dof.h"

using namespace std;

Solution::Solution(const System& system){
  // Save some data
  this->sol  = &(system.getSol());
  this->dofM = &(system.getDofManager());
  this->fs   = &(system.getFunctionSpace());

  // Get Mesh
  this->mesh = &(fs->getSupport().getMesh());  

  // Init
  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  // Scalar or Vector ?
  fsType = fs->getType();

  switch(fsType){
  case 0:
  case 3:
    isScalar = true;
    break;

  case 1:
  case 2:
    isScalar = false;
    break;
  }

  // Interpolate
  interpolate(fs);
}

Solution::~Solution(void){
  if(isScalar)
    delete nodalScalarValue;

  else
    delete nodalVectorValue;
}

void Solution::write(const std::string name, Writer& writer) const{
  // Get Domain
  const std::vector<const MElement*>& element
    = fs->getSupport().getAll();

  // Set Writer
  writer.setDomain(element);

  if(isScalar)
    writer.setValues(*nodalScalarValue);

  else
    writer.setValues(*nodalVectorValue);

  // Write
  writer.write(name);
}

void Solution::interpolate(const FunctionSpace* fs){
  // Init
  const unsigned int nTotVertex       = mesh->getVertexNumber();
  const std::vector<GroupOfDof*>& god = dofM->getAllGroups();
  const unsigned int nGod             = god.size();

  vector<bool> isInterpolated(nTotVertex, false);
  vector<MVertex*> node;


  // Scalar or Vector ?
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;

  if(isScalar){
    nodalScalarValue = new vector<double>(nTotVertex);
    fsScalar = static_cast<const FunctionSpaceScalar*>(fs);
  }

  else{
    nodalVectorValue = new vector<fullVector<double> >(nTotVertex);
    fsVector = static_cast<const FunctionSpaceVector*>(fs);
  }
    

  // Iterate on GroupOfElement
  for(unsigned int i = 0; i < nGod; i++){
    // Get Element
    MElement& element = 
      const_cast<MElement&>(god[i]->getGeoElement());

    // Get NodeS
    element.getVertices(node);
    const unsigned int nNode = node.size();

    // Iterate on Node
    for(unsigned int j = 0; j < nNode; j++){
      // Get *GMSH* Id
      const unsigned int id = node[j]->getNum() - 1;
      
      // If not interpolated: interpolate :-P !
      if(!isInterpolated[id]){
	// Get Dof
	const vector<const Dof*>& dof  = god[i]->getAll();
	const unsigned int        size = dof.size();

	// Get Coef
	vector<double> coef(size);
	for(unsigned int k = 0; k < size; k++){
	  // Look in Solution
	  coef[k] = 
	    (*sol)(dofM->getGlobalId(*dof[k])); 
	
	  // If Edge Space: Set Orientation
	  if(fsType == 1)
	    coef[k] *= god[i]->getOrientation(k);
	}
	
	// Get Node coordinate
	fullVector<double> xyz(3);
	xyz(0) = node[j]->x();
	xyz(1) = node[j]->y();
	xyz(2) = node[j]->z();

	// Interpolate (AT LAST !!)
	if(isScalar)
	  (*nodalScalarValue)[id] = 
	    fsScalar->interpolate(element, coef, xyz);

	else
	  (*nodalVectorValue)[id] = 
	    fsVector->interpolate(element, coef, xyz);

	isInterpolated[id] = true;
      }
    }
  }
}
