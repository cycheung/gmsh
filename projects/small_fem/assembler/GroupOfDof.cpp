#include "GroupOfDof.h"
#include <sstream>

GroupOfDof::GroupOfDof(const MElement& geoElement,
                       const std::vector<const Dof*>& dof,
                       const std::vector<size_t>& order){
  // Geo Element //
  element = &geoElement;

  // Alloc and Copy //
  this->dof   = new std::vector<const Dof*>(dof);
  this->order = new std::vector<size_t>(order);
}

GroupOfDof::GroupOfDof(const MElement& geoElement,
                       const std::vector<const Dof*>& dof){
  // Geo Element //
  element = &geoElement;

  // Alloc and Copy //
  this->dof = new std::vector<const Dof*>(dof);

  // Order //
  const size_t nDof = dof.size();

  this->order = new std::vector<size_t>(nDof);

  for(size_t i = 0; i < nDof; i++)
    (*this->order)[i] = i;
}

GroupOfDof::~GroupOfDof(void){
  delete dof;
  delete order;
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
