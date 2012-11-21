#include "GeoExtractor.h"

using namespace std;

map<const MVertex*, unsigned int, VertexComparator>*

GeoExtractor::extractVertex(const map<const MElement*, 
				      unsigned int, 
				      ElementComparator>& element){
  // Init //
  map<const MVertex*, unsigned int, VertexComparator>* 
    vertex = new map<const MVertex*, unsigned int, VertexComparator>;
  
  // Get Vertices //
  const map<const MElement*, unsigned int, ElementComparator>::const_iterator
    endE = element.end();
  
  map<const MElement*, unsigned int, ElementComparator>::const_iterator
    itE = element.begin();
  
  // Iterate on Elements
  for(; itE != endE; itE++){   
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Vertices
    const unsigned int N = myElement->getNumVertices();

    for(unsigned int i = 0; i < N; i++){
      // Take Current Vertex
      MVertex* myVertex = myElement->getVertex(i);

      // Try to Insert
      vertex->insert(pair<const MVertex* ,int>(myVertex, 0));
    }
  }

  // Return //
  return vertex;
}
