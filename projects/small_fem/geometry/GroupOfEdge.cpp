#include <set>
#include <sstream>
#include "GroupOfEdge.h"

using namespace std;

unsigned int GroupOfEdge::nextId = 0;

GroupOfEdge::GroupOfEdge(const GroupOfElement& goe){
  // Set Id
  id = nextId;
  nextId++;

  /*
  // Init Set //
  set<MEdge*, MEdgeLessThanNum> e;

  // Get Elements //
  const vector<MElement*>&  element = goe.getAll();  
  unsigned int             nElement = goe.getNumber();

  // Get Vertices //
  for(unsigned int i = 0; i < nElement; i++){
    vector<MEdge*> eEdge;
    element[i]->getVertices(eEdge);
    const unsigned int nEEdge = eEdge.size();

    for(unsigned int j = 0; j < nEEdge; j++)
      e.insert(eEdge[j]);
  }

  // Save Data //
  edge  = new vector<MEdge*>(e.begin(), e.end());
  nEdge = e.size();
  */
}

GroupOfEdge::~GroupOfEdge(void){
  delete edge;
}

string GroupOfEdge::toString(void) const{
  stringstream stream;
    
  stream << "*********************************************"    
	 << endl
	 << "* Group Of Edge #" << id   
	 << endl
	 << "*********************************************" 
	 << endl;
    
  return stream.str();
}
