#include "EdgeExtractor.h"

using namespace std;

EdgeExtractor::EdgeExtractor(void){
}


EdgeExtractor::~EdgeExtractor(void){
}

map<const MEdge*, unsigned int, EdgeComparator>*
EdgeExtractor::extract(const map<const MElement*, 
				 unsigned int, 
				 ElementComparator>& element){
  // Init //
  map<const MEdge*, unsigned int, EdgeComparator>* 
    edge = new map<const MEdge*, unsigned int, EdgeComparator>;
  
  // Get Edges //
  const map<const MElement*, unsigned int, ElementComparator>::const_iterator
    endE = element.end();
  
  map<const MElement*, unsigned int, ElementComparator>::const_iterator
    itE = element.begin();
  
  // Iterate on Elements
  for(; itE != endE; itE++){   
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Edges
    const unsigned int N = myElement->getNumEdges();

    for(unsigned int i = 0; i < N; i++){
      // Take Current Edge
      const MEdge myEdge = myElement->getEdge(i);
      
      // Make a copy (on heap)
      MEdge* edgeCopy = copy(myEdge);

      // Try to Insert
      pair<map<const MEdge*, unsigned int, EdgeComparator>::iterator,
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

MEdge* EdgeExtractor::copy(const MEdge& edge){
  MVertex* begin = edge.getVertex(0);
  MVertex* end   = edge.getVertex(1);

  return new MEdge(begin, end);  
}
