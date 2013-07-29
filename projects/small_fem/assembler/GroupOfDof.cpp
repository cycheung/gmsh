#include "GroupOfDof.h"
#include <sstream>

GroupOfDof::GroupOfDof(const MElement& element,
                       const std::vector<Dof>& dof){
  // Geo Element //
  this->element = &element;

  // Alloc and Copy //
  this->dof.assign(dof.begin(), dof.end());
}

GroupOfDof::~GroupOfDof(void){
}

std::string GroupOfDof::toString(void) const{
  std::stringstream stream;

  stream << "*************************** " << std::endl
         << "* Group Of Dof"               << std::endl
         << "*************************** " << std::endl
         << "* Associated Dofs:  "         << std::endl;

  for(size_t i = 0; i < dof.size(); i++)
    stream << "*    -- " << dof[i].toString() << std::endl;

  stream << "*************************** " << std::endl;

  return stream.str();
}
