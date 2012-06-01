#include <sstream>
#include "Dof.h"

Dof::Dof(const int entity, const int type){
  this->entity = entity;
  this->type = type;
}

std::string Dof::toString(void) const{
  std::stringstream s;
  
  s << "(" << entity << ", " << type << ")";

  return s.str();
}
