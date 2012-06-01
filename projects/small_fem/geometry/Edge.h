#ifndef _EDGE_H_
#define _EDGE_H_

#include <string>
#include "Entity.h"
#include "Node.h"

/**
   @class Edge
   @brief An Edge represents an edge of a Mesh
   
   Edge%s represent @em edges of a Mesh.@n

   An Edge is an Entity with 2 @em Nodes,
   namely '@c 0' and '@c 1'.

   @note
   As an Entity, an Edge can't be instantiated by
   the users
*/

class Mesh;

class Edge: public Entity{
 private:
  Node* node0;
  Node* node1;

  friend class Mesh;

 private:
  Edge(const int id, Node& node0, Node& node1);
  virtual ~Edge(void);

 public:
  Node& getNode0(void) const;
  Node& getNode1(void) const;

  virtual std::string toString(void) const;
};

/**
   @fn Edge::getNode0
   @return Returns the Node '@c 0' of this Edge

   @fn Edge::getNode1
   @return Returns the Node '@c 1' of this Edge
*/

//////////////////////
// Inline Functions //
//////////////////////

inline Edge::~Edge(void){
  // Edge is not responsible for deleting nodes
}

inline Node& Edge::getNode0(void) const{
  return *node0;
}

inline Node& Edge::getNode1(void) const{
  return *node1;
}

#endif
