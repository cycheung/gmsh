#include <set>
#include <sstream>
#include "GroupOfVertex.h"

using namespace std;

GroupOfVertex::GroupOfVertex(const GroupOfElement& goe){
  // Save Entity //
  this->id     = 0;

  vertex = new vector<MVertex*>(1);
}

GroupOfVertex::~GroupOfVertex(void){
  delete vertex;
}

string GroupOfVertex::toString(void) const{
  stringstream stream;
  
  return stream.str();
}
