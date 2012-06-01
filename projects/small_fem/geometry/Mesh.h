#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <string>
#include "Node.h"
#include "Edge.h"
#include "Element.h"

/**
   @class Mesh
   @brief Represents a mesh
   
   This class represents a mesh.@n

   A Mesh is composed of Entity, and of Element%s build 
   on these Entity.@n
 
   A Mesh is instantiated thanks to a 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   @c .msh file, wich discribes the mesh.@n

   @note
   Note that the Mesh will instantiate all its Entity,
   and its Element%s
*/

class Mesh{
 private:
  struct edgeTriplet{
    int   origin;
    int   end;
    Edge* edge;
  };
  
  class EdgeTripletComparator{
  public:
    inline bool operator()(const edgeTriplet& a, 
			   const edgeTriplet& b) const;
  };

  class PairComparator{
  public:
    inline bool operator()(const std::pair<int, int>& a, 
			   const std::pair<int, int>& b) const;
  };

  std::vector<Node*>* node;
  std::vector<Element*>* nodeElement; 
  std::vector<Element*>* edgeElement;
  std::vector<Edge*>* edge;
  int nNode;
  int nNodeElement;
  int nLine;
  int nEdge;
  int nEdgeElement;

  std::ifstream* in;  
  
 public:
   Mesh(const std::string fileName);
  ~Mesh(void);
   
  Node& getNode(const int i) const;
  const std::vector<Node*>& getAllNodes(void) const;
  int getNbNode(void) const;

  Element& getNodeElement(const int i) const;
  const std::vector<Element*>& getAllNodeElements(void) const;
  int getNbNodeElement(void) const;

  Edge& getEdge(const int i) const;
  const std::vector<Edge*>& getAllEdges(void) const;
  int getNbEdge(void) const;

  Element& getEdgeElement(const int i) const;
  const std::vector<Element*>& getAllEdgeElements(void) const;
  int getNbEdgeElement(void) const;

  std::string toString(void) const;

 private:
  void skipLine(const int N);
  bool parse(const std::string str);

  void getTri(int id, std::istringstream& stream);
  void getLine(int id, std::istringstream& stream);
  void getEdges(void);
};

/**
   @fn Mesh::Mesh
   @param fileName The path to the @c .msh file discribing the Meh
   @return Returns a new Mesh

   @fn Mesh::~Mesh
   @return Deletes this Mesh

   @fn Mesh::getNode
   @param i The @c ID of a Node in the Mesh
   @return Returns the requested Node

   @fn Mesh::getAllNodes
   @return Returns the set of all Node%s in the Mesh

   @fn Mesh::getNbNode
   @return Returns the number of Node%s in the Mesh

   @fn Mesh::getNodeElement
   @param i The @c ID of a @em nodal Element
   in the Mesh
   @return Returns the requested @em nodal Element

   @fn Mesh::getAllNodeElements
   @return Returns the set of all 
   @em nodal Element%s in the Mesh

   @fn Mesh::getNbNodeElement
   @return Returns the number of 
   @em nodal Element%s in the Mesh

   @fn Mesh::getEdge
   @param i The @c ID of an Edge in the Mesh
   @return Returns the requested Edge

   @fn Mesh::getAllEdges
   @return Returns the set of all Edge%s in the Mesh

   @fn Mesh::getNbEdge
   @return Returns the number of Edge%s in the Mesh

   @fn Mesh::getEdgeElement
   @param i The @c ID of an @em edge Element
   in the Mesh
   @return Returns the requested @em edge Element

   @fn Mesh::getAllEdgeElements
   @return Returns the set of all 
   @em edge Element%s in the Mesh

   @fn Mesh::getNbEdgeElement
   @return Returns the number of 
   @em edge Element%s in the Mesh

   @fn Mesh::toString
   @return Returns a string associated to the Mesh
 */


//////////////////////
// Inline Functions //
//////////////////////

inline Node& Mesh::getNode(const int i) const{
  return *((*node)[i]);
}
  
inline const std::vector<Node*>& Mesh::getAllNodes(void) const{
  return *node;
}

inline int Mesh::getNbNode(void) const{
  return nNode;
}

inline Element& Mesh::getNodeElement(const int i) const{
  return *((*nodeElement)[i]);
}
  
inline const std::vector<Element*>& Mesh::getAllNodeElements(void) const{
  return *nodeElement;
}

inline int Mesh::getNbNodeElement(void) const{
  return nNodeElement; 
}

inline Edge& Mesh::getEdge(const int i) const{
  return *((*edge)[i]);
}

inline const std::vector<Edge*>& Mesh::getAllEdges(void) const{
  return *edge;
}

inline int Mesh::getNbEdge(void) const{
  return nEdge;
}

inline Element& Mesh::getEdgeElement(const int i) const{
  return *((*edgeElement)[i]);
}

inline const std::vector<Element*>& Mesh::getAllEdgeElements(void) const{
  return *edgeElement;
}

inline int Mesh::getNbEdgeElement(void) const{
  return nEdgeElement;
}

inline bool Mesh::PairComparator::operator()(const std::pair<int, int>& a, 
				       const std::pair<int, int>& b) const{
  return (a.first < b.first) || 
    ((a.first == b.first) && (a.second < b.second));
}

inline bool Mesh::EdgeTripletComparator::operator()(const edgeTriplet& a, 
					  const edgeTriplet& b) const{
  return (a.origin < b.origin) || 
    ((a.origin == b.origin) && (a.end < b.end));
}


#endif
