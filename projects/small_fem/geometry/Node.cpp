#include <sstream>
#include "Node.h"

Node::Node(const int id, const double x, const double y, const double z):
  Entity::Entity(id, 0){
  
  this->x  = x;
  this->y  = y;
  this->z  = z;
}

std::string Node::toString(void) const{
  std::stringstream stream;
  
  stream << "Node " << id << std::endl;
  stream << "**** " << std::endl;
  stream << "  -- X: " << x << std::endl;
  stream << "  -- Y: " << y << std::endl;
  stream << "  -- Z: " << z << std::endl;
  if(hasPhysical)
    stream << "  -- Physical: " << physical << std::endl;
  if(hasValue)
    stream << "  -- Value: " << value << std::endl;

  return stream.str();
}
