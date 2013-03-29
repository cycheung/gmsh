#include <sstream>
#include "Dof.h"

Dof::Dof(void){
  this->entity  = 0;
  this->type    = 0;
}

Dof::Dof(const Dof& other){
  this->entity  = other.entity;
  this->type    = other.type;
}

Dof::Dof(unsigned long int entity, unsigned int type){
  this->entity  = entity;
  this->type    = type;
}

Dof::~Dof(void){
}

std::string Dof::toString(void) const{
  std::stringstream s;

  s << "(" << entity << ", " << type << ")";

  return s.str();
}
