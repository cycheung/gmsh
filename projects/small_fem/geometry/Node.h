#ifndef _NODE_H_
#define _NODE_H_

#include <string>
#include "Entity.h"

/**
   @class Node
   @brief A Node represents a node of a Mesh
   
   Node%s represent @em verticies of a Mesh.@n

   A Node is an Entity with 3 @em coordinates,
   namely @c X, @c Y and @c Z.

   @note
   As an Entity, a Node can't be instantiated by
   the users
*/

class Mesh;

class Node: public Entity{
 private:
  double x;
  double y;
  double z;

  friend class Mesh;

 private:
  Node(const int id, const double x, const double y, const double z);
  virtual ~Node(void);

 public:
  double getX(void) const;
  double getY(void) const;
  double getZ(void) const;

  virtual std::string toString(void) const;
};

/**
   @fn Node::getX
   @return Returns the @c X coordinate of this Node

   @fn Node::getY
   @return Returns the @c Y coordinate of this Node

   @fn Node::getZ
   @return Returns the @c Z coordinate of this Node

*/

//////////////////////
// Inline Functions //
//////////////////////

inline Node::~Node(void){
}

inline double Node::getX(void) const{
  return x;
}

inline double Node::getY(void) const{
  return y;
}

inline double Node::getZ(void) const{
  return z;
}

#endif
