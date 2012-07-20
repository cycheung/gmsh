#include <sstream>
#include <set>

#include "Polynomial.h"
#include "PlotBasis.h"

using namespace std;

PlotBasis::PlotBasis(const GroupOfElement& group, const Basis& basis){
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
      delete (*nodalScalarValue)[i];
   
    delete nodalScalarValue;
  }

  if(nodalVectorValue){
    for(int i = 0; i < nFunction; i++)
      delete (*nodalVectorValue)[i];
    
    delete nodalVectorValue;
  }
}

void PlotBasis::write(const string name) const{
  for(int i = 0; i < nFunction; i++){
    stringstream nameCat, fileName; 
    
    nameCat  << name << i + 1;
    fileName << nameCat.str() << ".msh";

    ofstream out;
    out.open(fileName.str().c_str());
    
    writeHeader(out);
    writeNodes(out);
    writeElements(out);
    
    writeNodalValues(out, nameCat.str(), i);  

    out.close();
  }
}

void PlotBasis::writeHeader(ofstream& out) const{
  out << "$MeshFormat" << endl;
  out << "2.2 0 8" << endl;
  out << "$EndMeshFormat" << endl; 
}

void PlotBasis::writeNodes(ofstream& out) const{
  out << "$Nodes" << endl;
  out << N << endl;

  for(int i = 0; i < N; i++){
    out << (*node)[i]->getNum() << " ";
    out << (*node)[i]->x()      << " ";
    out << (*node)[i]->y()      << " ";
    out	<< (*node)[i]->z()      << endl;
  }

  out << "$EndNodes" << endl;
}

void PlotBasis::writeElements(ofstream& out) const{
  out << "$Elements" << endl;
  out << E << endl;
  
  for(int i = 0; i < E; i++){
    out << (*element)[i]->getNum()        << " " 
	<< (*element)[i]->getTypeForMSH() 
	<< " 2 1 1" << " "; 
           // 2 Tags --> (1 physical entity, 1 elementary geometry) 

    const int M = (*element)[i]->getNumVertices();
    for(int j = 0; j < M; j++)
      out << (*element)[i]->getVertex(j)->getNum() << " ";
    out << endl;
  }

  out << "$EndElements" << endl;
}

void PlotBasis::writeNodalValues(ofstream& out, 
				 const string name,
				 int fun) const{

  out << "$ElementNodeData"   << endl
      << "1"                  << endl  // 1 string tag
      << "\"" << name << "\"" << endl  // (name)
      << "1"                  << endl  // 1 real tag 
      << "0"                  << endl  // (time value)
      << "3"                  << endl  // 3 integer tag
      << "0"                  << endl; // (time step index)
    
  if(isScalar)
    out << "1" << endl;                // (number of field -- scalar)
  else
    out << "3" << endl;                // (number of field -- vector)
    
  out << E << endl;                    // (number of element)
  
  for(int i = 0; i < E; i++){
    out << (*element)[i]->getNum()         << " " 
	<< (*element)[i]->getNumVertices() << " ";
    
    const int M = (*element)[i]->getNumVertices();
    for(int j = 0; j < M; j++){
      const int id = (*element)[i]->getVertex(j)->getNum() - 1;
      // Note: getNum() ranges from *1* to MAX
      //   --> we need to substract 1 !!

      if(isScalar)
	out << (*(*nodalScalarValue)[fun])[id] << " ";
      else
	out << (*(*nodalVectorValue)[fun])[id](0) << " "
	    << (*(*nodalVectorValue)[fun])[id](1) << " "
	    << (*(*nodalVectorValue)[fun])[id](2) << " ";
    }

    out << endl;
  }
  out << "$EndElementNodeData" << endl;
}

void PlotBasis::getGeometry(const GroupOfElement& group){
  // Get All Elements //
  element = &(group.getAll());
  E       = element->size();

  // Get All Vertices //
  set<MVertex*, MVertexLessThanNum> setVertex;

  for(int i = 0; i < E; i++){
    const int N = (*element)[i]->getNumVertices();
    
    for(int j = 0; j < N; j++)
      setVertex.insert((*element)[i]->getVertex(j));
  }

  // Serialize the set into a vector //
  node = new vector<MVertex*>(setVertex.begin(), 
			      setVertex.end());
  N    = node->size();
}


void PlotBasis::interpolate(const BasisScalar& basis){
  // Allocate //
  nodalVectorValue = NULL;  
  nodalScalarValue = new vector<vector<double>*>(nFunction);

  for(int i = 0; i < nFunction; i++){
    (*nodalScalarValue)[i] = new vector<double>(N);
  }

  // Get Functions //
  const vector<Polynomial>& fun = basis.getBasis();
  
  // Interpolate //
  for(int f = 0; f < nFunction; f++){
    for(int n = 0; n < N; n++){
      (*(*nodalScalarValue)[f])[n] = 
	fun[f].at((*node)[n]->x(),
		  (*node)[n]->y(),
		  (*node)[n]->z());
    }
  }  
}

void PlotBasis::interpolate(const BasisVector& basis){
  // Allocate //
  nodalScalarValue = NULL;  
  nodalVectorValue = new vector<vector<fullVector<double> >*>(nFunction);

  for(int i = 0; i < nFunction; i++){
    (*nodalVectorValue)[i] = new vector<fullVector<double> >(N);
  }

  // Get Functions //
  const vector<vector<Polynomial> >& fun = basis.getBasis();
  
  // Interpolate //
  for(int f = 0; f < nFunction; f++){
    for(int n = 0; n < N; n++){
      (*(*nodalVectorValue)[f])[n] = 
	Polynomial::at(fun[f], 
		       (*node)[n]->x(),
		       (*node)[n]->y(),
		       (*node)[n]->z());
    }
  }
}
