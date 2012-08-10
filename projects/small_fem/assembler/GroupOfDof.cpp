#include "GroupOfDof.h"

#include <sstream>

unsigned int GroupOfDof::nextId = 0;

GroupOfDof::GroupOfDof(unsigned int numberOfDof, const MElement& geoElement,
		       const FunctionSpace& fs , const Mesh& mesh){
  // Set ID //
  id = nextId;
  nextId++;

  // Get FS, Mesh and Element //
  this->element = &geoElement;
  this->fs      = &fs;
  this->mesh    = &mesh;

  // Set GroupOfDof //
  nDof = numberOfDof;
  dof  = new std::vector<const Dof*>(nDof);

  nextDof = 0;
}

GroupOfDof::~GroupOfDof(void){
  // GroupOfDofs are not responsible for
  // deleting dofs, orientations and Jacobian
  delete dof;
}

int GroupOfDof::getOrientation(unsigned int dofId) const{
  // Take Requested Dof //
  const Dof& dof = *(*(this->dof))[dofId];

  // Get Type // 
  // Return 0 If *NOT* an Edge //
  if(fs->getElementType(dof) != 1)
    return 0;

  // If we have an Edge //
  // Get *Oriented* Edge
  const MEdge& oEdge = mesh->getEdge(fs->getElementGlobalId(dof));

  // And Look for its orientation in the *Geometric* Element
  return orientation(*element, oEdge);
}

void GroupOfDof::add(const Dof* dof){
  this->dof->at(nextDof) = dof;
  nextDof++;
}

int GroupOfDof::orientation(const MElement& element, 
			    const MEdge& edge){
  
  MElement& eelement = 
    const_cast<MElement&>(element);

  const unsigned int nEdge = 
    eelement.getNumEdges();

  for(unsigned int i = 0; i < nEdge; i++)
    if(equal(edge, eelement.getEdge(i)))
      return 1;

  return -1;
}

std::string GroupOfDof::toString(void) const{
  std::stringstream stream;
  
  stream << "*************************** " << std::endl
	 << "* Group Of Dof #" << id       << std::endl
	 << "*************************** " << std::endl
	 << "* Associated Dofs:  " << std::endl;

  for(unsigned int i = 0; i < nDof; i++)
    stream << "*    -- " << get(i).toString() << std::endl;

  stream << "*************************** " << std::endl;

  return stream.str();
}
