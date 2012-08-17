#include <sstream>
#include <set>

#include "BasisScalar.h"
#include "BasisVector.h"
#include "Polynomial.h"

#include "PlotBasis.h"

using namespace std;

PlotBasis::PlotBasis(const Basis& basis,
		     const GroupOfElement& group, 
		     Writer& writer){

  this->writer = &writer;
  
  nFunction = basis.getSize();
  getGeometry(group);

  if(basis.isScalar()){
    isScalar = true;
    interpolate(static_cast<const BasisScalar&>(basis));
  }

  else{
    isScalar = false;
    interpolate(static_cast<const BasisVector&>(basis));
  }
}

PlotBasis::~PlotBasis(void){
  delete node;

  if(nodalScalarValue){
    for(int i = 0; i < nFunction; i++)
      delete nodalScalarValue[i];
   
    delete[] nodalScalarValue;
  }

  if(nodalVectorValue){
    for(int i = 0; i < nFunction; i++)
      delete nodalVectorValue[i];
    
    delete[] nodalVectorValue;
  }
}

void PlotBasis::plot(const string name) const{
  for(int i = 0; i < nFunction; i++){
    stringstream nameCat; 
    nameCat  << name << i + 1;
    
    if(isScalar)
      writer->setValues(*(nodalScalarValue[i]));

    else
      writer->setValues(*(nodalVectorValue[i]));

    writer->write(nameCat.str());
  }
}

void PlotBasis::getGeometry(const GroupOfElement& group){
  // Get All Elements //
  element = &(group.getAll());
  E       = element->size();

  // Get All Vertices //
  set<MVertex*, MVertexLessThanNum> setVertex;

  for(int i = 0; i < E; i++){
    const int N = (*element)[i]->getNumVertices();
    MElement* myElement = 
      const_cast<MElement*>((*element)[i]);
    
    for(int j = 0; j < N; j++)
      setVertex.insert(myElement->getVertex(j));
  }

  // Serialize the set into a vector //
  node = new vector<MVertex*>(setVertex.begin(), 
			      setVertex.end());
  N    = node->size();
}


void PlotBasis::interpolate(const BasisScalar& basis){
  // Allocate //
  nodalVectorValue = NULL;
  nodalScalarValue = new vector<double>*[nFunction];

  for(int i = 0; i < nFunction; i++)
    nodalScalarValue[i] = new vector<double>(N);

  // Get Functions //
  const vector<Polynomial>& fun = basis.getFunctions();
  
  // Interpolate //
  for(int f = 0; f < nFunction; f++){
    for(int n = 0; n < N; n++){
      (*nodalScalarValue[f])[n] = 
	fun[f].at((*node)[n]->x(),
		  (*node)[n]->y(),
		  (*node)[n]->z());
    }
  }  
}

void PlotBasis::interpolate(const BasisVector& basis){
  // Allocate //
  nodalScalarValue = NULL;  
  nodalVectorValue = new vector<fullVector<double> >*[nFunction];

  for(int i = 0; i < nFunction; i++)
    nodalVectorValue[i] = new vector<fullVector<double> >(N);

  // Get Functions //
  const vector<vector<Polynomial> >& fun = basis.getFunctions();
  
  // Interpolate //
  for(int f = 0; f < nFunction; f++){
    for(int n = 0; n < N; n++){
      (*nodalVectorValue[f])[n] = 
	Polynomial::at(fun[f], 
		       (*node)[n]->x(),
		       (*node)[n]->y(),
		       (*node)[n]->z());
    }
  }
}
