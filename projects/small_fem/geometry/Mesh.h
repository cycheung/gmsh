#ifndef _MESH_H_
#define _MESH_H_

#include <map>
#include <string>

#include "Comparators.h"
#include "GModel.h"

#include "GroupOfElement.h"
#include "MElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   This class is responsible of the handling mesh elements
   (Such as Quads, Tets, Edges, Vertices, ...).@n

   It is also responsible of the @em numbering of those
   elements.@n
   Indeed, each element shall be granted a @em unique @c ID.@n

   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   .@c msh file, wich discribes the mesh.@n
*/

class GroupOfElement;

class Mesh{
 private:
  GModel* model;

  std::map<const MElement*, unsigned int, ElementComparator>* element;           
  std::map<const MVertex*, unsigned int, VertexComparator>*   vertex;
  std::map<const MEdge*, unsigned int, EdgeComparator>*       edge;
  std::map<const MFace*, unsigned int, FaceComparator>*       face;  

  std::map<unsigned int, const MElement*>* idElement;
  std::map<unsigned int, const MVertex*>*  idVertex;
  std::map<unsigned int, const MEdge*>*    idEdge;
  std::map<unsigned int, const MFace*>*    idFace;

  std::multimap<int, const MElement*>* physical;

  int nextId;

 public:
   Mesh(const std::string fileName);
  ~Mesh(void);

  GModel& getModel(void) const;

  unsigned int getGlobalId(const MElement& element) const;
  unsigned int getGlobalId(const MVertex& vertex) const;
  unsigned int getGlobalId(const MEdge& edge) const;
  unsigned int getGlobalId(const MFace& face) const;

  const MElement& getElement(unsigned int id) const;
  const MVertex&  getVertex(unsigned int id) const;
  const MEdge&    getEdge(unsigned int id) const;
  const MFace&    getFace(unsigned int id) const;

  const std::vector<const MVertex*> getAllVertex(void) const;

  unsigned int getElementNumber(void) const;
  unsigned int getVertexNumber(void) const;
  unsigned int getEdgeNumber(void) const;
  unsigned int getFaceNumber(void) const;

  GroupOfElement getFromPhysical(int physicalId) const;
  
  std::string toString(void) const;
 
 private:
  void number(void);
};


/**
   @fn Mesh::Mesh
   @param fileName The path to the @c .msh file discribing the Mesh

   Instanciates a new Mesh
   **
   
   @fn Mesh::~Mesh
   Deletes this Mesh
   **

   @fn Mesh::getModel
   @return Returns the Model used for generating this Mesh
   **

   @fn unsigned int Mesh::getGlobalId(const MElement& element) const
   @param element A MElement
   @return Returns the @em global @em @c ID (in this Mesh) of the
   given MElement
   **

   @fn unsigned int Mesh::getGlobalId(const MVertex& vertex) const
   @param vertex A MVertex
   @return Returns the @em global @em @c ID (in this Mesh) of the
   given MVertex
   **

   @fn unsigned int Mesh::getGlobalId(const MEdge& edge) const
   @param edge A MEdge
   @return Returns the @em global @em @c ID (in this Mesh) of the
   given MEdge
   **

   @fn unsigned int Mesh::getGlobalId(const MFace& face) const
   @param face A MFace
   @return Returns the @em global @em @c ID (in this Mesh) of the
   given MFace
   **
 
   @fn Mesh::getElement
   @param id A natural number
   @return Returns the MElement with the given @em global @c ID
   @note 
   If no MElement has the given global @c ID, an Exception is thrown
   **

   @fn Mesh::getVertex
   @param id A natural number
   @return Returns the MVertex with the given @em global @c ID
   @note 
   If no MVertex has the given global @c ID, an Exception is thrown
   **

   @fn Mesh::getEdge
   @param id A natural number
   @return Returns the MEdge with the given @em global @c ID
   @note 
   If no MEdge has the given global @c ID, an Exception is thrown
   **

   @fn Mesh::getFace
   @param id A natural number
   @return Returns the MFace with the given @em global @c ID
   @note 
   If no MFace has the given global @c ID, an Exception is thrown
   **

   @fn Mesh::getAllVertex
   @return Returns all the Vertices of this Mesh
   **
   
   @fn Mesh::getElementNumber
   @return Returns the number of Element in this Mesh
   @note By Element we mean Quads, Tets, etc@n
   This excludes Vertices, Edges, Faces and Cells
   **

   @fn Mesh::getVertexNumber
   @return Returns the number of MVertices in this Mesh
   **

   @fn Mesh::getEdgeNumber
   @return Returns the number of MEdge%s in this Mesh
   **

   @fn Mesh::getFaceNumber
   @return Returns the number of MFace%s in this Mesh
   **

   @fn Mesh::getFromPhysical
   @param physicalId A physical @c ID 
   (see <a href="http://www.geuz.org/gmsh">gmsh</a> 
   documentation)

   @return @em Instantiate a new GroupOfElement, containing
   the MElements of the given physical @c ID
   **

   @fn Mesh::toString
   @return Returns a description of this Mesh
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline GModel& Mesh::getModel(void) const{
  return *model;
}

inline unsigned int Mesh::getElementNumber(void) const{
  return element->size();
}

inline unsigned int Mesh::getVertexNumber(void) const{
  return vertex->size();
}

inline unsigned int Mesh::getEdgeNumber(void) const{
  return edge->size();
}

inline unsigned int Mesh::getFaceNumber(void) const{
  return 0;
}

#endif
