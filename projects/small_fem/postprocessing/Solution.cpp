#include "Solution.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "Exception.h"
#include "GModel.h"
#include "MElement.h"
#include "MVertex.h"
#include "Dof.h"

using namespace std;

void Solution::initSystem(const System& system){
  // Save some data
  this->dofM = &(system.getDofManager());
  this->fs   = &(system.getFunctionSpace());

  // Get Mesh
  this->mesh = &(fs->getSupport().getMesh());  

  // Get Solution 
  ownSol = false;
  sol    = &(system.getSol());

  // Init
  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  // Scalar or Vector ?
  fsType = fs->getType();

  switch(fsType){
  case 0:
  case 3:
    scalar = true;
    break;

  case 1:
  case 2:
    scalar = false;
    break;
  }
}

void Solution::initSystem(const EigenSystem& system, 
			  unsigned int eigenNumber){
  // Save some data
  this->dofM = &(system.getDofManager());
  this->fs   = &(system.getFunctionSpace());

  // Get Mesh
  this->mesh = &(fs->getSupport().getMesh());  

  // Get Solution
  ownSol = true;
  sol    = getSol(system.getEigenVectors(), eigenNumber);
    
  // Init
  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  // Scalar or Vector ?
  fsType = fs->getType();

  switch(fsType){
  case 0:
  case 3:
    scalar = true;
    break;

  case 1:
  case 2:
    scalar = false;
    break;
  }
}

Solution::Solution(const System& system){
  // Init
  initSystem(system);

  // Get Visu Domain
  this->visuDomain = &(fs->getSupport());

  // Interpolate 
  // NB: interpolate() is faster than 
  // interpolateOnVisu() (no Octree)
  interpolate();
}

Solution::Solution(const System& system,
		   const GroupOfElement& visu){
  // Init
  initSystem(system);

  // Get Visu Domain
  this->visuDomain = &visu;

  // Interpolate
  // NB: Can't use interpolate(), because
  // we don't have *all* nodes (visu mesh) 
  // with Dofs
  interpolateOnVisu();
}

Solution::Solution(const EigenSystem& system,
		   unsigned int eigenNumber){
  // Init
  initSystem(system, eigenNumber);

  // Get Visu Domain
  this->visuDomain = &(fs->getSupport());

  // Interpolate 
  // NB: interpolate() is faster than 
  // interpolateOnVisu() (no Octree)
  interpolate();
}

Solution::Solution(const EigenSystem& system,
		   unsigned int eigenNumber,
		   const GroupOfElement& visu){
  // Init
  initSystem(system, eigenNumber);

  // Get Visu Domain
  this->visuDomain = &visu;

  // Interpolate
  // NB: Can't use interpolate(), because
  // we don't have *all* nodes (visu mesh) 
  // with Dofs
  interpolateOnVisu();
}

Solution::Solution(double (*f)(fullVector<double>& xyz), 
		   const GroupOfElement& visu){
  // Init
  scalar = true;
  ownSol = false;
  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  // Get Visu Domain
  this->visuDomain = &visu;

  // Get Function
  fScalar = f;

  // Evaluate f on visu
  evaluateF();
}

Solution::Solution(fullVector<double> (*f)(fullVector<double>& xyz), 
		   const GroupOfElement& visu){
  // Init
  scalar = false;
  ownSol = false;
  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  // Get Visu Domain
  this->visuDomain = &visu;

  // Get Function
  fVector = f;

  // Evaluate f on visu
  evaluateF();
}

Solution::~Solution(void){
  if(ownSol)
    delete sol;

  if(scalar)
    delete nodalScalarValue;

  else
    delete nodalVectorValue;
}

void Solution::write(const std::string name, Writer& writer) const{
  // Set Writer
  writer.setDomain(visuDomain->getAll());

  if(scalar)
    writer.setValues(*nodalScalarValue);

  else
    writer.setValues(*nodalVectorValue);

  // Write
  writer.write(name);
}

void Solution::interpolate(void){
  // Init
  const unsigned int nTotVertex       = mesh->getVertexNumber();
  const std::vector<GroupOfDof*>& god = fs->getAllGroups();
  const unsigned int nGod             = god.size();

  vector<bool> isInterpolated(nTotVertex, false);
  vector<MVertex*> node;


  // Scalar or Vector ?
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;

  if(scalar){
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
	for(unsigned int k = 0; k < size; k++)
	  // Look in Solution
	  coef[k] = 
	    (*sol)(dofM->getGlobalId(*dof[k])); 

	// Get Node coordinate
	fullVector<double> xyz(3);
	xyz(0) = node[j]->x();
	xyz(1) = node[j]->y();
	xyz(2) = node[j]->z();

	// Interpolate (AT LAST !!)
	if(scalar)
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

void Solution::interpolateOnVisu(void){
  // Init
  const Mesh& visuMesh              = visuDomain->getMesh();
  const unsigned int nTotVertex     = visuMesh.getVertexNumber();
  const vector<const MVertex*> node = visuMesh.getAllVertex(); 

  // Scalar or Vector ?
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;

  if(scalar){
    nodalScalarValue = new vector<double>(nTotVertex);
    fsScalar = static_cast<const FunctionSpaceScalar*>(fs);
  }

  else{
    nodalVectorValue = new vector<fullVector<double> >(nTotVertex);
    fsVector = static_cast<const FunctionSpaceVector*>(fs);
  }
    
  // Get Model for Octrees
  GModel& model = mesh->getModel();
  const int dim = model.getDim();

  // Iterate on *NODES*
  for(unsigned int i = 0; i < nTotVertex; i++){
    // Search element (in System Mesh) containg this 
    // visu node
    SPoint3   point   = node[i]->point();
    MElement* element = model.getMeshElementByCoord(point, dim, true);
    
    // WARNING: if no element found, set to zero
    if(!element){
      if(scalar)
	(*nodalScalarValue)[node[i]->getNum() - 1] = 0;
    
      else
	(*nodalVectorValue)[node[i]->getNum() - 1] = 0;
    }

    else{
      // Get GroupOfDof related to this Element
      const GroupOfDof& god = fs->getGoDFromElement(*element);
      
      // Get Dof
      const vector<const Dof*>& dof  = god.getAll();
      const unsigned int        size = dof.size();
      
      // Get Coef
      vector<double> coef(size);
      for(unsigned int k = 0; k < size; k++)
	// Look in Solution
	coef[k] = 
	  (*sol)(dofM->getGlobalId(*dof[k])); 
      
      // Get Node coordinate
      fullVector<double> xyz(3);
      xyz(0) = node[i]->x();
      xyz(1) = node[i]->y();
      xyz(2) = node[i]->z();
      
      // Interpolate (AT LAST !!)
      if(scalar)
	(*nodalScalarValue)[node[i]->getNum() - 1] = 
	  fsScalar->interpolate(*element, coef, xyz);
      
      else
	(*nodalVectorValue)[node[i]->getNum() - 1] = 
	  fsVector->interpolate(*element, coef, xyz);
    }
  }
}

void Solution::evaluateF(void){
  // Init
  const Mesh& visuMesh              = visuDomain->getMesh();
  const unsigned int nTotVertex     = visuMesh.getVertexNumber();
  const vector<const MVertex*> node = visuMesh.getAllVertex(); 

  // Scalar or Vector ?
  if(scalar)
    nodalScalarValue = new vector<double>(nTotVertex);

  else
    nodalVectorValue = new vector<fullVector<double> >(nTotVertex);
    
  // Iterate on *NODES*
  for(unsigned int i = 0; i < nTotVertex; i++){
    // Get Node coordinate
    fullVector<double> xyz(3);
    xyz(0) = node[i]->x();
    xyz(1) = node[i]->y();
    xyz(2) = node[i]->z();

    // Evaluate (AT LAST !!)
    if(scalar)
      (*nodalScalarValue)[node[i]->getNum() - 1] = 
	fScalar(xyz);
    
    else
      (*nodalVectorValue)[node[i]->getNum() - 1] = 
	fVector(xyz);
  }
}

const fullVector<double>* Solution::
getSol(const vector<vector<complex<double> > >& eVector,
       unsigned int eigenNumber){
  
  // Init 
  unsigned int size       = eVector[eigenNumber].size(); 
  fullVector<double>* sol = new fullVector<double>(size);

  // Get Sol
  for(unsigned int i = 0; i < size; i++){
    (*sol)(i) = real(eVector[eigenNumber][i]);
  }

  // Return
  return sol;
}

vector<double>& Solution::getNodalScalarValue(void) const{
  if(!scalar)
    throw Exception("Solution: try to get Scalar value in a Vectorial Solution");

  else
    return *nodalScalarValue;
}

vector<fullVector<double> >& Solution::getNodalVectorValue(void) const{
  if(scalar)
    throw Exception("Solution: try to get Vectorial value in a Scalar Solution");

  else
    return *nodalVectorValue;  
}
