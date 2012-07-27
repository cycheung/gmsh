#include <sstream>

#include "Dof.h"
#include "Exception.h"

Dof::Dof(const unsigned int entity, const unsigned int type){
  this->entity  = entity;
  this->type    = type;
  this->unknown = true;
}

Dof::~Dof(void){
}

double Dof::getValue(void) const{
  if(unknown)
    throw Exception
      ("Dof %s is an unknown, and has no value", 
       toString().c_str());

  return value;
}

std::string Dof::toString(void) const{
  std::stringstream s;
  
  s << "(" << entity << ", " << type << ")";

  return s.str();
}
