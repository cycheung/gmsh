#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <string>

#include "Group.h"
#include "GEntity.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   A Mesh is composed of @em physicals and of @em Group%s @em of @em elements @n
   defined on those physicals.

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n
*/

class Mesh{
 private:
  
 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  int    getNbPhysicals(void) const;
  Group& getGroup(int i) const;

  const std::vector<int>&    getAllPhysicals(void) const;
  const std::vector<Group*>& getAllGroups(void) const;

  std::string toString(void) const;
};

/**
   @fn Mesh::Mesh
   @param fileName The path to the @c .msh file discribing the Mesh
   @return Returns a new Mesh

   @fn Mesh::~Mesh
   @return Deletes this Mesh
 */


//////////////////////
// Inline Functions //
//////////////////////


#endif
