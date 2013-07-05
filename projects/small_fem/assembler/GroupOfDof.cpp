#include "GroupOfDof.h"
#include <sstream>

GroupOfDof::GroupOfDof(const MElement& geoElement,
                       const std::vector<const Dof*>& dof){
  // Geo Element //
  element = &geoElement;

  // Alloc and Copy //
  this->dof = new std::vector<const Dof*>(dof);
}

GroupOfDof::~GroupOfDof(void){
  delete dof;
}

std::string GroupOfDof::toString(void) const{
  std::stringstream stream;

  stream << "*************************** " << std::endl
         << "* Group Of Dof"               << std::endl
         << "*************************** " << std::endl
         << "* Associated Dofs:  "         << std::endl;

  for(size_t i = 0; i < dof->size(); i++)
    stream << "*    -- " << (*dof)[i]->toString() << std::endl;

  stream << "*************************** " << std::endl;

  return stream.str();
}
