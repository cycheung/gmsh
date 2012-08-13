#include <set>
#include "Writer.h"

using namespace std;

Writer::Writer(void){
  hasValue  = false;
  hasDomain = false;
  node      = NULL;
}

Writer::~Writer(void){
  if(node)
    delete node;
}

void Writer::setValues(const std::vector<double>& value){
  hasValue = true;
  isScalar = true;

  nodalScalarValue = &value;
  nodalVectorValue = NULL;
}

void Writer::setValues(const std::vector<fullVector<double> >& value){
  hasValue = true;
  isScalar = false;

  nodalScalarValue = NULL;
  nodalVectorValue = &value;
}

void Writer::setDomain(const std::vector<const MElement*>& element){
  // Erease Old Domain (if one) //
  if(hasDomain)
    delete node;

  // Get Elements //
  this->element = &element;
  this->E       = element.size();
  
  // Get All Vertices //
  set<MVertex*, MVertexLessThanNum> setVertex;

  for(int i = 0; i < E; i++){
    const int N = element[i]->getNumVertices();
    MElement* myElement = 
      const_cast<MElement*>(element[i]);

    for(int j = 0; j < N; j++)
      setVertex.insert(myElement->getVertex(j));
  }

  // Serialize the set into a vector //
  node = new vector<MVertex*>(setVertex.begin(), 
			      setVertex.end());
  N    = node->size();

  // Set hasDomain //
  hasDomain = true;
}
