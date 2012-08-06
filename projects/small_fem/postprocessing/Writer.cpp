#include <set>
#include "Writer.h"

using namespace std;

Writer::Writer(void){
  hasValue  = false;
  hasDomain = false;
}

Writer::~Writer(void){
  delete node;
}

void Writer::setValues(std::vector<double>& value){
  hasValue = true;
  isScalar = true;

  nodalScalarValue = &value;
  nodalVectorValue = NULL;
}

void Writer::setValues(std::vector<fullVector<double> >& value){
  hasValue = true;
  isScalar = false;

  nodalScalarValue = NULL;
  nodalVectorValue = &value;
}

void Writer::setDomain(const std::vector<MElement*>& element){
  // Get Elements //
  this->element = &element;
  this->E       = element.size();
  
  // Set hasDomain
  hasDomain = true;
  
  // Get All Vertices //
  set<MVertex*, MVertexLessThanNum> setVertex;

  for(int i = 0; i < E; i++){
    const int N = element[i]->getNumVertices();
    
    for(int j = 0; j < N; j++)
      setVertex.insert(element[i]->getVertex(j));
  }

  // Serialize the set into a vector //
  node = new vector<MVertex*>(setVertex.begin(), 
			      setVertex.end());
  N    = node->size();
}
