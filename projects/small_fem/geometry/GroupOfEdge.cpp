#include <set>
#include <list>
#include <sstream>

#include "GroupOfEdge.h"
#include "Exception.h"

using namespace std;

unsigned int GroupOfEdge::nextId = 0;

GroupOfEdge::GroupOfEdge(const GroupOfElement& goe, Mesh& mesh){
  // Set Id
  id = nextId;
  nextId++;

  // Init //
  this->mesh = &mesh;
  orientation = new map<MEdge, int, EdgeComparator>;
  
  // Get Elements //
  const vector<MElement*>&  element = goe.getAll();  
  unsigned int             nElement = goe.getNumber();

  // Get Edges //
  for(unsigned int i = 0; i < nElement; i++){
    const unsigned int N = element[i]->getNumEdges();

    for(unsigned int j = 0; j < N; j++){
      // Take Current Edge
      MEdge myEdge = element[i]->getEdge(j);
      
      // Try to Insert
      pair<map<MEdge, int, EdgeComparator>::iterator,
	   bool> insert = 
	orientation->insert(pair<MEdge ,int>(myEdge, 1));

      // If Insertion is a success,
      // Insert inverted Edge
      if(insert.second)
	orientation->insert
	  (pair<MEdge ,int>(invert(myEdge), -1));	
    }
  }

  // Keep Edge With good Orientation //
  map<MEdge, int, EdgeComparator>::iterator end = 
    orientation->end();

  map<MEdge, int, EdgeComparator>::iterator it = 
    orientation->begin();

  unsigned int i = 0;

  list<MEdge*> lst;

  for(; it != end; i++, it++)
    if(it->second == 1)
      lst.push_back(const_cast<MEdge*>(&(it->first)));

  // Edges Vector //
  edge  = new vector<MEdge*>(lst.begin(), lst.end());
  nEdge = edge->size(); 
}

GroupOfEdge::~GroupOfEdge(void){
  delete orientation;
  delete edge;
}

int GroupOfEdge::getOrientation(const MEdge& edge) const{
  // Iterator //
  const map<MEdge, int, EdgeComparator>::iterator end = 
    orientation->end();
  
  map<MEdge, int, EdgeComparator>::iterator it = 
    orientation->begin();

  // Check If we have the requested Edge
  if(it != end)
    return it->second;
  
  // Else: Exception !
  else
    throw Exception("Can't find Edge for Orientation");
}

string GroupOfEdge::toString(void) const{
  stringstream stream;

  const map<MEdge, int, EdgeComparator>::iterator end = 
    orientation->end();
  
  map<MEdge, int, EdgeComparator>::iterator it = 
    orientation->begin();
    
  stream << "*********************************************"    
	 << endl
	 << "* Group Of Edge #" << id   
	 << endl
	 << "*********************************************" 
	 << endl<< "*" 
	 << endl
    	 << "* This group contains the following edges: " << endl;


  for(; it != end; it++)
    stream << "*    -- ID: " 
	   << mesh->getGlobalId(it->first) 
	   << " (Orientation: " << it->second << ")" 
	   << endl;


  stream << "*********************************************" 
	 << endl;

  return stream.str();
}

MEdge GroupOfEdge::invert(MEdge& edge){
  MVertex* begin = edge.getVertex(0);
  MVertex* end   = edge.getVertex(1);

  return MEdge(end, begin);
}
