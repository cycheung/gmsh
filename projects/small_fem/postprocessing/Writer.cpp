#include <set>
#include "Writer.h"

using namespace std;

Writer::Writer(void){
  ownSol    = false;
  hasValue  = false;
  hasDomain = false;
  node      = NULL;
}

Writer::~Writer(void){
  if(ownSol)
    delete sol;

  if(node)
    delete node;
}

void Writer::setValues(const std::vector<double>& value){
  nodalScalarValue = &value;
  nodalVectorValue = NULL;

  fs   = NULL;
  dofM = NULL;
  sol  = NULL;

  ownSol   = false;
  hasValue = true;
  isScalar = true;
  isNodal  = true;
}

void Writer::setValues(const std::vector<fullVector<double> >& value){
  nodalScalarValue = NULL;
  nodalVectorValue = &value;

  fs   = NULL;
  dofM = NULL;
  sol  = NULL;

  ownSol   = false;
  hasValue = true;
  isScalar = false;
  isNodal  = true;
}

void Writer::setValues(const System& system){
  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  fs   = &(system.getFormulation().fs());
  dofM = &(system.getDofManager());
  sol  = &(system.getSol());

  if(system.isSolved()){
    ownSol   = false;
    hasValue = true;
    isScalar = fs->isScalar();
    isNodal  = false;
  }

  else{
    ownSol   = false;
    hasValue = false;
  }

  setDomain(fs->getSupport().getAll());
}

void Writer::setValues(const EigenSystem& system, unsigned int eigenNumber){
  // Delete old Sol if any, and if Own
  if(ownSol)
    delete sol;

  nodalScalarValue = NULL;
  nodalVectorValue = NULL;

  fs   = &(system.getFormulation().fs());
  dofM = &(system.getDofManager());
  sol  = getSol(system.getEigenVectors(), eigenNumber);

  if(system.isSolved()){    
    ownSol   = true;
    hasValue = true;
    isScalar = fs->isScalar();
    isNodal  = false;
  }

  else{
    hasValue = false;
  }

  setDomain(fs->getSupport().getAll());
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

const fullVector<double>* Writer::
getSol(const vector<vector<complex<double> > >& eVector,
       unsigned int eigenNumber){
  
  // Init 
  unsigned int size       = eVector[eigenNumber].size(); 
  fullVector<double>* sol = new fullVector<double>(size);

  // Get Sol
  for(unsigned int i = 0; i < size; i++){
    (*sol)(i) = real(eVector[eigenNumber][i]);
  }

  // Return
  return sol;
}
