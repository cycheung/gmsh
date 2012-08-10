#include "Solution.h"
#include "MVertex.h"
#include "Dof.h"

using namespace std;

Solution::Solution(const System& system, const Mesh& mesh){
  // Get Mesh
  this->mesh = &mesh;

  // Get FunctionSpace
  fs = &(system.getFunctionSpace());

  // Save some data
  this->sol     = &(system.getSol());
  this->dofM    = &(system.getDofManager());
  this->god     = &(dofM->getAllGroups());
  this->element = &(fs->getSupport().getAll());

  // Interpolate
  switch(fs->getType()){
  case 0:
  case 3:
    interpolateScalar(static_cast<const FunctionSpaceScalar&>(*fs));
    break;

  case 1:
  case 2:
    interpolateVector(static_cast<const FunctionSpaceVector&>(*fs));
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

void Solution::interpolateScalar(const FunctionSpaceScalar& fs){
  isScalar = true;

  const vector<const Dof*> dof = dofM->getAllDofs();
  unsigned int nDof            = dof.size();

  nodalScalarValue = new vector<double>(nDof);

  for(unsigned int i = 0; i < nDof; i++){
    (*nodalScalarValue)
      [mesh->getVertex
       (fs.getElementGlobalId(*dof[i])).getNum() - 1] = 
      
      (*sol)(dofM->getGlobalId(*dof[i]));
  }
}

void Solution::interpolateVector(const FunctionSpaceVector& fs){
  // Init
  isScalar = false;
  const unsigned int nDof = dofM->dofNumber();
  const unsigned int nGod = god->size();

  nodalVectorValue = new vector<fullVector<double> >(nDof);
  vector<bool> isInterpolated(nDof, false);
  vector<MVertex*> node;
  
  // Iterate on GroupOfElement
  for(unsigned int i = 0; i < nGod; i++){
    // Get Element
    MElement& element = 
      const_cast<MElement&>((*god)[i]->getGeoElement());

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
	const vector<const Dof*>& dof  = (*god)[i]->getAll();
	const unsigned int        size = dof.size();

	// Get Coef
	vector<double> coef(size);
	for(unsigned int k = 0; k < size; k++)
	  coef[k] = 
	    (*sol)(dofM->getGlobalId(*dof[k])) * 
	    (*god)[i]->getOrientation(k);
	

	// Get Node coordinate
	fullVector<double> xyz(3);
	xyz(0) = node[j]->x();
	xyz(1) = node[j]->y();
	xyz(2) = node[j]->z();

	// Interpolate (AT LAST !!)
	(*nodalVectorValue)[id] = fs.interpolate(element, coef, xyz);
	isInterpolated[id]      = true;
      }
    }
  }
}
