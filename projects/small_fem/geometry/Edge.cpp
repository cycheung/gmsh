#include <sstream>
#include "Edge.h"

Edge::Edge(const int id, Node& node0, Node& node1):
  Entity::Entity(id, 1){

  this->node0 = &node0;
  this->node1 = &node1;  
}

std::string Edge::toString(void) const{
  std::stringstream stream;
  
  stream << "Edge " << id << std::endl;
  stream << "**** " << std::endl;
  stream << "  -- Node0: " << node0->getId() << std::endl;
  stream << "  -- Node1: " << node1->getId() << std::endl;
  if(hasPhysical)
    stream << "  -- Physical: " << physical << std::endl;
  if(hasValue)
    stream << "  -- Value: " << value << std::endl;

  return stream.str();
}
