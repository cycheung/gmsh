#include <sstream>
#include "Entity.h"

Entity::Entity(const int id, const int type){
  this->id   = id;
  this->type = type;
  
  hasValue    = false;
  hasPhysical = false;
}

std::string Entity::toString(void) const{
  std::stringstream stream;
  
  stream << "Entity " << id << std::endl;
  stream << "****** " << std::endl;

  if(hasPhysical)
    stream << "  -- Physical: " << physical << std::endl;
  if(hasValue)
    stream << "  -- Value: " << value << std::endl;

  return stream.str();
}
