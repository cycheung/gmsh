#include <set>
#include <sstream>
#include "GroupOfVertex.h"

using namespace std;

GroupOfVertex::GroupOfVertex(const GroupOfElement& goe){
  // Init Set //
  set<MVertex*, MVertexLessThanNum> v;

  // Get Elements //
  const vector<MElement*>&  element = goe.getAll();  
  unsigned int             nElement = goe.getNumber();

  // Get Vertices //
  for(unsigned int i = 0; i < nElement; i++){
    vector<MVertex*> eVertex;
    element[i]->getVertices(eVertex);
    const unsigned int nEVertex = eVertex.size();

    for(unsigned int j = 0; j < nEVertex; j++)
      v.insert(eVertex[j]);
  }

  // Save Data //
  vertex  = new vector<MVertex*>(v.begin(), v.end());
  nVertex = v.size();
}

GroupOfVertex::~GroupOfVertex(void){
  delete vertex;
}

string GroupOfVertex::toString(void) const{
  stringstream stream;
  
  stream << "*********************************************"    
	 << endl
	 << "* Group Of Vertex                          *"    
	 << endl
	 << "*********************************************" 
	 << endl << "*" 
	 << endl
	 << "* This group contains the following vertices: " << endl;

  for(unsigned int i = 0; i < nVertex; i++){
    stream << "*    -- ID: " 
	   << (*vertex)[i]->getNum() << endl;
  }
  
  stream << "*********************************************" 
	 << endl;
  
  return stream.str();
}
