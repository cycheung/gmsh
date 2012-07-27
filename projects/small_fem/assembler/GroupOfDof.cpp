#include <sstream>
#include "GroupOfDof.h"

GroupOfDof::GroupOfDof(int numberOfDof, const MElement& geoElement){
  element = &geoElement;
  
  nDof = numberOfDof;
  dof  = new std::vector<Dof*>(nDof);

  nextDof = 0;
}

GroupOfDof::~GroupOfDof(void){
  // GroupOfDofs are not responsible for
  // deleting dofs, orientations and Jacobian
  delete dof;
}

void GroupOfDof::add(Dof* dof){
  this->dof->at(nextDof) = dof;
  nextDof++;
}

std::string GroupOfDof::toString(void) const{
  std::stringstream stream;
  
  stream << "******************* " << std::endl
	 << "* GroupOfDof number " << getId() << std::endl
	 << "******************* " << std::endl
	 << "* Associated Dofs:  " << std::endl;

  for(int i = 0; i < nDof; i++){
    stream << "*    -- " << get(i).toString() << std::endl;
    stream << "*      Entity: " << const_cast<MElement*>(element)->getVertex(i)->getNum()
	   << std::endl;
  }

  stream << "******************* " << std::endl;

  return stream.str();
}
