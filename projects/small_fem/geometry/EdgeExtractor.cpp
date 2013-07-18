#include "GeoExtractor.h"

using namespace std;

map<const MEdge*, size_t, EdgeComparator>*
GeoExtractor::extractEdge(const map<const MElement*,
                                    size_t,
                          ElementComparator>& element){
  // Init //
  map<const MEdge*, size_t, EdgeComparator>*
    edge = new map<const MEdge*, size_t, EdgeComparator>;

  // Get Edges //
  const map<const MElement*, size_t, ElementComparator>::const_iterator
    endE = element.end();

  map<const MElement*, size_t, ElementComparator>::const_iterator
    itE = element.begin();

  // Iterate on Elements
  for(; itE != endE; itE++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Edges
    const size_t N = myElement->getNumEdges();

    for(size_t i = 0; i < N; i++){
      // Take Current Edge
      const MEdge myEdge = myElement->getEdge(i);

      // Make a copy (on heap)
      MEdge* edgeCopy = copy(myEdge);

      // Try to Insert
      pair<map<const MEdge*, size_t, EdgeComparator>::iterator,
           bool> insert =
        edge->insert(pair<const MEdge* ,int>(edgeCopy, 0));

      // If Insertion is not a success,
      // Delete edgeCopy
      if(!insert.second)
        delete edgeCopy;
    }
  }

  // Return //
  return edge;
}
