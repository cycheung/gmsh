#include <vector>
#include <sstream>

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
    map<const MElement*, unsigned int, ElementComparator>*,
    multimap<int, const MElement*>*
    >
    elementsExtracted = ElementExtractor::extract(entity);

  element  = elementsExtracted.first;
  physical = elementsExtracted.second;

  // Extract Nodes //
  vertex = VertexExtractor::extract(*element);

  // Extract Edges //
  pair<
    map<const MEdge*, unsigned int, EdgeComparator>*, 
    map<const MEdge*, int, OrientedEdgeComparator>*
    > 
    edgesExtracted = EdgeExtractor::extract(*element);
    
  edge        = edgesExtracted.first;
  orientation = edgesExtracted.second;

  // Number Mesh Entities //
  idVertex  = new map<unsigned int, const MVertex*>;
  idEdge    = new map<unsigned int, const MEdge*>;
  idElement = new map<unsigned int, const MElement*>;

  nextId = 0;
  number();
}

Mesh::~Mesh(void){
  // Delete Elements //
  
  // WARNING
  // Mesh is *NOT* responsible for
  // Deleting MElement*
  delete idElement;
  delete element;
  delete physical;

  // Delete Vertices

  // WARNING
  // Mesh is *NOT* responsible for
  // Deleting MVertex*
  delete idVertex;
  delete vertex;

  // Delete Edges //
  const map<const MEdge*, int, OrientedEdgeComparator>::iterator
    endE = orientation->end();
  
  map<const MEdge*, int, OrientedEdgeComparator>::iterator 
    itE = orientation->begin();

  for(; itE != endE; itE++)
    delete itE->first;

  delete idEdge;
  delete edge;  
  delete orientation;

  // Delete Model //
  delete model;
}

unsigned int Mesh::getGlobalId(const MElement& element) const{
  map<const MElement*, 
      unsigned int, 
      ElementComparator>::iterator it = 
    this->element->find(&element);

  if(it == this->element->end())
    throw 
      Exception("Element not found");

  return it->second;   
}

unsigned int Mesh::getGlobalId(const MVertex& vertex) const{
  map<const MVertex*, 
      unsigned int, 
      VertexComparator>::iterator it = 
    this->vertex->find(&vertex);

  if(it == this->vertex->end())
    throw 
      Exception("Vertex not found");

  return it->second; 
}

unsigned int Mesh::getGlobalId(const MEdge& edge) const{
  // WARNING:
  // Here, we use a Edge Comparator,
  // such that Edge Orientation
  // do *NOT* matter !!

  // Look for Edge //
  map<const MEdge*, 
      unsigned int, 
      EdgeComparator>::iterator it = 
    this->edge->find(&edge);

  if(it == this->edge->end()){
    throw 
      Exception("Edge not found");
  }

  return it->second; 
}

const MElement& Mesh::getElement(unsigned int id) const{
  map<unsigned int, const MElement*>::iterator it = 
    idElement->find(id);

  if(it == idElement->end())
    throw 
      Exception("No Element with Global Id %d found", id);

  return *(it->second);   
}

const MVertex& Mesh::getVertex(unsigned int id) const{
  map<unsigned int, const MVertex*>::iterator it = 
    idVertex->find(id);

  if(it == idVertex->end())
    throw 
      Exception("No Vertex with Global Id %d found", id);

  return *(it->second);   
}

const MEdge& Mesh::getEdge(unsigned int id) const{
  map<unsigned int, const MEdge*>::iterator it = 
    idEdge->find(id);

  if(it == idEdge->end())
    throw 
      Exception("No Edge with Global Id %d found", id);

  return *(it->second);   
}

void Mesh::number(void){
  // Get Iterators //
  const map<const MElement*, unsigned int, ElementComparator>::iterator
    endEl = element->end();           

  const map<const MVertex*, unsigned int, VertexComparator>::iterator
    endV = vertex->end();           

  const map<const MEdge*, unsigned int, ElementComparator>::iterator
    endEd = edge->end();           
  
  map<const MElement*, unsigned int, ElementComparator>::iterator
    itEl = element->begin();

  map<const MVertex*, unsigned int, VertexComparator>::iterator
    itV = vertex->begin();           
  
  map<const MEdge*, unsigned int, ElementComparator>::iterator
    itEd = edge->begin();
  
  // Number Vertices //
  for(; itV != endV; itV++){
    itV->second = nextId;
    idVertex->insert
      (pair<unsigned int, const MVertex*>
       (itV->second, itV->first));
    
    nextId++;
  }

  // Number Edges //
  for(; itEd != endEd; itEd++){
    itEd->second = nextId;
    idEdge->insert
      (pair<unsigned int, const MEdge*>
       (itEd->second, itEd->first));
    
    nextId++;
  }

  // Number Elements //
  for(; itEl != endEl; itEl++){
    itEl->second = nextId;
    idElement->insert
      (pair<unsigned int, const MElement*>
       (itEl->second, itEl->first));
    
    nextId++;
  }
}

int Mesh::getOrientation(const MEdge& edge) const{
  // Look for Edge //
  map<const MEdge*, 
      int, 
      OrientedEdgeComparator>::iterator it = 
    orientation->find(&edge);

  if(it == orientation->end()){
    throw 
      Exception("Edge not found");
  }

  return it->second; 
}

GroupOfElement Mesh::getFromPhysical(int physicalId) const{
  const std::pair<std::multimap<int, const MElement*>::iterator, 
                  std::multimap<int, const MElement*>::iterator> p = 
    physical->equal_range(physicalId);
  
  return GroupOfElement(p.first, p.second, *this);
}

string Mesh::toString(void) const{
  // Iterators //
  const map<const MElement*, unsigned int, ElementComparator>::iterator
    endEl = element->end();           

  const map<const MVertex*, unsigned int, VertexComparator>::iterator
    endV = vertex->end();           

  const map<const MEdge*, unsigned int, ElementComparator>::iterator
    endEd = edge->end();           
  
  map<const MElement*, unsigned int, ElementComparator>::iterator
    itEl = element->begin();

  map<const MVertex*, unsigned int, VertexComparator>::iterator
    itV = vertex->begin();           
  
  map<const MEdge*, unsigned int, ElementComparator>::iterator
    itEd = edge->begin();
  
  stringstream stream;
  

  // Header //
  stream << "***********************************************"    
	 << endl
	 << "*                     Mesh                    *"    
	 << endl
	 << "***********************************************"
	 << endl; 


  // Elements //
  stream << "*                                             *"
	 << endl
	 << "* This mesh contains the following Elements:  *" 
	 << endl;
  
  for(; itEl != endEl; itEl++)
    stream << "*   -- Element "
	   << getGlobalId(*itEl->first)
	   << endl;

  stream << "*                                             *"
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
	   << getGlobalId(*itV->first)
	   << endl
	   << "*    (["
	   << itV->first->x()
	   << ", "
	   << itV->first->y()
	   << ", "
	   << itV->first->z()
	   << "])"
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
	   << endl
	   << "*    (["
	   << itEd->first->getVertex(0)->x()
	   << ", "
	   << itEd->first->getVertex(0)->y()
	   << ", "
	   << itEd->first->getVertex(0)->z()
	   << "], ["
	   << itEd->first->getVertex(1)->x()
	   << ", "
	   << itEd->first->getVertex(1)->y()
	   << ", "
	   << itEd->first->getVertex(1)->z()
	   << "])"
	   << endl;


  stream << "*                                             *"
	 << endl
	 << "***********************************************"  
	 << endl;

  
  // Retrun //
  return stream.str();
}
