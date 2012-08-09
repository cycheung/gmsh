#include "EdgeExtractor.h"

using namespace std;

EdgeExtractor::EdgeExtractor(void){
}


EdgeExtractor::~EdgeExtractor(void){
}

pair<map<const MEdge*, unsigned int, EdgeComparator>*, 
     map<const MEdge*, int, OrientedEdgeComparator>*> 

EdgeExtractor::extract(const map<const MElement*, 
				 unsigned int, 
				 ElementComparator>& element){
  // Init //
  map<const MEdge*, int, OrientedEdgeComparator>* 
    orientation = new map<const MEdge*, int, OrientedEdgeComparator>;

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
      pair<map<const MEdge*, int, OrientedEdgeComparator>::iterator,
	   bool> insert = 
	orientation->insert(pair<const MEdge* ,int>(edgeCopy, 1));

      // If Insertion is a success,
      // Insert inverted Edge
      if(insert.second)
	orientation->insert
	  (pair<const MEdge* ,int>(invert(myEdge), -1));

      // Else, Delete edgeCopy
      else
	delete edgeCopy;
    }
  }

  // Keep Edge With orientation of +1 //
  map<const MEdge*, int, OrientedEdgeComparator>::iterator end = 
    orientation->end();

  map<const MEdge*, int, OrientedEdgeComparator>::iterator it = 
    orientation->begin();

  for(; it != end; it++)
    if(it->second == 1)
      edge->insert
	(pair<const MEdge*, unsigned int>(it->first, 0));

  // Return //
  return 
    pair<map<const MEdge*, unsigned int, EdgeComparator>*, 
	 map<const MEdge*, int, OrientedEdgeComparator>*> 
  (edge, orientation);
}

MEdge* EdgeExtractor::copy(const MEdge& edge){
  MVertex* begin = edge.getVertex(0);
  MVertex* end   = edge.getVertex(1);

  return new MEdge(begin, end);  
}

MEdge* EdgeExtractor::invert(const MEdge& edge){
  MVertex* begin = edge.getVertex(0);
  MVertex* end   = edge.getVertex(1);

  return new MEdge(end, begin);
}
