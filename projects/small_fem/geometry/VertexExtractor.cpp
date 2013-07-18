#include "GeoExtractor.h"

using namespace std;

map<const MVertex*, size_t, VertexComparator>*
GeoExtractor::extractVertex(const map<const MElement*,
                                      size_t,
                                      ElementComparator>& element){
  // Init //
  map<const MVertex*, size_t, VertexComparator>*
    vertex = new map<const MVertex*, size_t, VertexComparator>;

  // Get Vertices //
  const map<const MElement*, size_t, ElementComparator>::const_iterator
    endE = element.end();

  map<const MElement*, size_t, ElementComparator>::const_iterator
    itE = element.begin();

  // Iterate on Elements
  for(; itE != endE; itE++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Vertices
    const size_t N = myElement->getNumVertices();

    for(size_t i = 0; i < N; i++){
      // Take Current Vertex
      MVertex* myVertex = myElement->getVertex(i);

      // Try to Insert
      vertex->insert(pair<const MVertex* ,int>(myVertex, 0));
    }
  }

  // Return //
  return vertex;
}
