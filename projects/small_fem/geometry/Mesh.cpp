#include <vector>
#include <sstream>

#include "GeoExtractor.h"
#include "Exception.h"
#include "Mesh.h"

using namespace std;

Mesh::Mesh(const std::string fileName){
  // New Mode //
  model = new GModel("SmallFEM");

  // Read Mesh //
  if(!model->readMSH(fileName))
    throw Exception("Can't open file: %s", fileName.c_str());

  // Get Entity //
  vector<GEntity*> entity;
  model->getEntities(entity);

  // Extract Element //
  pair<
    map<const MElement*, size_t, ElementComparator>*,
    multimap<int, const MElement*>*
    >
    elementsExtracted = GeoExtractor::extractElement(entity);

  element  = elementsExtracted.first;
  physical = elementsExtracted.second;

  // Extract Nodes, Edges and Faces //
  vertex = GeoExtractor::extractVertex(*element);
  edge   = GeoExtractor::extractEdge(*element);
  face   = GeoExtractor::extractFace(*element);

  // Number Geometry //
  nextId = 0;
  number();
}

Mesh::~Mesh(void){
  // Delete Elements //

  // WARNING
  // Mesh is *NOT* responsible for
  // Deleting MElement*
  delete element;
  delete physical;

  // Delete Vertices //

  // WARNING
  // Mesh is *NOT* responsible for
  // Deleting MVertex*
  delete vertex;

  // Delete Edges //
  const map<const MEdge*, size_t, EdgeComparator>::iterator
    endE = edge->end();

  map<const MEdge*, size_t, EdgeComparator>::iterator
    itE = edge->begin();

  for(; itE != endE; itE++)
    delete itE->first;

  delete edge;

  // Delete Faces //
  const map<const MFace*, size_t, FaceComparator>::iterator
    endF = face->end();

  map<const MFace*, size_t, FaceComparator>::iterator
    itF = face->begin();

  for(; itF != endF; itF++)
    delete itF->first;

  delete face;

  // Delete Model //
  delete model;
}

size_t Mesh::getGlobalId(const MElement& element) const{
  map<const MElement*,
      size_t,
      ElementComparator>::iterator it =
    this->element->find(&element);

  if(it == this->element->end())
    throw
      Exception("Element not found");

  return it->second;
}

size_t Mesh::getGlobalId(const MVertex& vertex) const{
  map<const MVertex*,
      size_t,
      VertexComparator>::iterator it =
    this->vertex->find(&vertex);

  if(it == this->vertex->end())
    throw
      Exception("Vertex not found");

  return it->second;
}

size_t Mesh::getGlobalId(const MEdge& edge) const{
  // Look for Edge //
  map<const MEdge*,
      size_t,
      EdgeComparator>::iterator it =
    this->edge->find(&edge);

  if(it == this->edge->end()){
    throw
      Exception("Edge not found");
  }

  return it->second;
}

size_t Mesh::getGlobalId(const MFace& face) const{
  // Look for Face //
  map<const MFace*,
      size_t,
      FaceComparator>::iterator it =
    this->face->find(&face);

  if(it == this->face->end()){
    throw
      Exception("Face not found");
  }

  return it->second;
}

const vector<const MVertex*> Mesh::getAllVertex(void) const{
  // Init
  const size_t size = vertex->size();

  map<const MVertex*, size_t, VertexComparator>::iterator
    itV = vertex->begin();

  // Alloc
  vector<const MVertex*> v(size);

  // Fill Vector
  for(size_t i = 0; i < size; i++, itV++)
    v[i] = itV->first;

  // Return
  return v;
}

void Mesh::getAllVertexCoordinate(fullMatrix<double>& coord) const{
  // Get All Vertex
  const vector<const MVertex*>  vertex = getAllVertex();
  const size_t                 nVertex = vertex.size();

  // Allocate 'coord' and populate
  coord.resize(nVertex, 3);

  for(size_t i = 0; i < nVertex; i++){
    coord(i, 0) = vertex[i]->x();
    coord(i, 1) = vertex[i]->y();
    coord(i, 2) = vertex[i]->z();
  }
}

void Mesh::number(void){
  // Get Iterators //
  const map<const MElement*, size_t, ElementComparator>::iterator
    endEl = element->end();

  const map<const MVertex*, size_t, VertexComparator>::iterator
    endV = vertex->end();

  const map<const MEdge*, size_t, EdgeComparator>::iterator
    endEd = edge->end();

  const map<const MFace*, size_t, FaceComparator>::iterator
    endF = face->end();

  map<const MElement*, size_t, ElementComparator>::iterator
    itEl = element->begin();

  map<const MVertex*, size_t, VertexComparator>::iterator
    itV = vertex->begin();

  map<const MEdge*, size_t, EdgeComparator>::iterator
    itEd = edge->begin();

  map<const MFace*, size_t, FaceComparator>::iterator
    itF = face->begin();

  // Number Vertices //
  for(; itV != endV; itV++){
    itV->second = nextId;
    nextId++;
  }

  // Number Edges //
  for(; itEd != endEd; itEd++){
    itEd->second = nextId;
    nextId++;
  }

  // Number Faces //
  for(; itF != endF; itF++){
    itF->second = nextId;
    nextId++;
  }

  // Number Elements //
  for(; itEl != endEl; itEl++){
    itEl->second = nextId;
    nextId++;
  }
}

GroupOfElement Mesh::getFromPhysical(int physicalId) const{
  const std::pair<std::multimap<int, const MElement*>::iterator,
                  std::multimap<int, const MElement*>::iterator> p =
    physical->equal_range(physicalId);

  return GroupOfElement(p.first, p.second, *this);
}

string Mesh::toString(void) const{
  // Iterators //
  const map<const MElement*, size_t, ElementComparator>::iterator
    endEl = element->end();

  const map<const MVertex*, size_t, VertexComparator>::iterator
    endV = vertex->end();

  const map<const MEdge*, size_t, EdgeComparator>::iterator
    endEd = edge->end();

  const map<const MFace*, size_t, FaceComparator>::iterator
    endF = face->end();

  map<const MElement*, size_t, ElementComparator>::iterator
    itEl = element->begin();

  map<const MVertex*, size_t, VertexComparator>::iterator
    itV = vertex->begin();

  map<const MEdge*, size_t, EdgeComparator>::iterator
    itEd = edge->begin();

  map<const MFace*, size_t, FaceComparator>::iterator
    itF = face->begin();

  stringstream stream;


  // Header //
  stream << "***********************************************"
         << endl
         << "*                     Mesh                    *"
         << endl
         << "***********************************************"
         << endl;


  // Vertices //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Vertex:    *"
         << endl;

  for(; itV != endV; itV++)
    stream << "*   -- Vertex "
           << showpos
           << getGlobalId(*itV->first)
           << ": {"
           << itV->first->x()
           << ", "
           << itV->first->y()
           << ", "
           << itV->first->z()
           << "}"
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;


  // Edges //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Edges:     *"
         << endl;

  for(; itEd != endEd; itEd++)
    stream << "*   -- Edge "
           << getGlobalId(*itEd->first)
           << ": ["
           << getGlobalId(*itEd->first->getVertex(0))
           << ", "
           << getGlobalId(*itEd->first->getVertex(1))
           << "]"
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  // Faces //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Faces:     *"
         << endl;

  for(; itF != endF; itF++)
    stream << "*   -- Face "
           << getGlobalId(*itF->first)
           << ": ["
           << getGlobalId(*itF->first->getVertex(0))
           << ", "
           << getGlobalId(*itF->first->getVertex(1))
           << ", "
           << getGlobalId(*itF->first->getVertex(2))
           << "]"
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;


  // Elements //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Elements:  *"
         << endl;

  for(; itEl != endEl; itEl++){
    int nVertex =
      const_cast<MElement*>(itEl->first)->getNumPrimaryVertices();

    int nVertexMinus = nVertex - 1;

    stream << "*   -- Element "
           << getGlobalId(*itEl->first)
           << ": [";

    for(int i = 0; i < nVertexMinus; i++)
      stream << getGlobalId(*const_cast<MElement*>(itEl->first)->getVertex(i))
             << ", ";

    stream <<
      getGlobalId(*const_cast<MElement*>(itEl->first)->getVertex(nVertexMinus))
           << "]"
           << endl;
  }

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  // Retrun //
  return stream.str();
}
