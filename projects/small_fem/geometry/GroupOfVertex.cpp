#include <set>
#include <sstream>
#include "GroupOfVertex.h"

using namespace std;

unsigned int GroupOfVertex::nextId = 0;

GroupOfVertex::GroupOfVertex(const GroupOfElement& goe, Mesh& mesh){
  // Set Id //
  id = nextId;
  nextId++;

  // Init //
  this->mesh = &mesh;
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
	 << "* Group Of Vertex #" << id   
	 << endl
	 << "*********************************************" 
	 << endl << "*" 
	 << endl
	 << "* This group contains the following vertices: " << endl;

  for(unsigned int i = 0; i < nVertex; i++)
    stream << "*    -- ID: " 
	   << mesh->getGlobalId(*(*vertex)[i]) << endl;
  
  stream << "*********************************************" 
	 << endl;
  
  return stream.str();
}
