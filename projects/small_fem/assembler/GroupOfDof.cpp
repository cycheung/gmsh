#include "GroupOfDof.h"

#include <sstream>

unsigned int GroupOfDof::nextId = 0;

GroupOfDof::GroupOfDof(unsigned int numberOfDof, const MElement& geoElement){
  // Init //
  element = &geoElement;
  id      = nextId;
  nextId++;

  // Set GroupOfDof //
  nDof = numberOfDof;
  dof  = new std::vector<const Dof*>(nDof);

  nextDof = 0;
}

GroupOfDof::~GroupOfDof(void){
  delete dof;
}

void GroupOfDof::add(const Dof* dof){
  this->dof->at(nextDof) = dof;
  nextDof++;
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
