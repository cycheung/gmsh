#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <map>
#include <string>

#include "Comparators.h"
#include "GroupOfElement.h"
#include "GModel.h"

#include "ElementExtractor.h"
#include "VertexExtractor.h"
#include "EdgeExtractor.h"

#include "MElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   A Mesh is a collection of @em GroupOfElement%s 

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n
*/

class Mesh{
 private:
  GModel* model;

  std::map<const MElement*, unsigned int, ElementComparator>* element;           
  std::map<const MVertex*, unsigned int, MVertexLessThanNum>* vertex;
  std::map<const MEdge*, unsigned int, EdgeComparator>*       edge;

  std::map<unsigned int, const MVertex*>*  idVertex;
  std::map<unsigned int, const MEdge*>*    idEdge;
  std::map<unsigned int, const MElement*>* idElement;           

  std::multimap<int, const MElement*>*                 physical;
  std::map<const MEdge*, int, OrientedEdgeComparator>* orientation;

  int nextId;

 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  unsigned int getGlobalId(const MElement& element) const;
  unsigned int getGlobalId(const MVertex& vertex) const;
  unsigned int getGlobalId(const MEdge& edge) const;
  unsigned int getGlobalId(const MFace& face) const;

  const MElement& getElement(unsigned int id) const;
  const MVertex&  getVertex(unsigned int id) const;
  const MEdge&    getEdge(unsigned int id) const;
  const MFace&    getFace(unsigned int id) const;

  GroupOfElement getFromPhysical(int physicalId) const;
  
  std::string toString(void) const;
 
 private:
  void number(void);
};


/**
   @fn Mesh::Mesh
   Instanciate a new Mesh
   @param fileName The path to the @c .msh file discribing the Mesh
   
   @fn Mesh::~Mesh
   Deletes this Mesh
 
   @fn Mesh::toString
   @return Returns a description of this Mesh
*/


//////////////////////
// Inline Functions //
//////////////////////

inline GroupOfElement Mesh::getFromPhysical(int physicalId) const{
  const std::pair<std::multimap<int, const MElement*>::iterator, 
                  std::multimap<int, const MElement*>::iterator> p = 
    physical->equal_range(physicalId);
  
  return GroupOfElement(p.first, p.second);
}

#endif
