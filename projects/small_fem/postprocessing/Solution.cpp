#include "Solution.h"
#include "Formulation.h"
#include "InterpolatorScalar.h"
#include "InterpolatorVector.h"

using namespace std;

Solution::Solution(const Mesh& mesh, const Formulation& formulation){
  /*
  msh     = &mesh;
  element = &msh->getAllNodeElements();
  E       = msh->getNbNodeElement();
  node    = &msh->getAllNodes();
  N       = msh->getNbNode();
  
  interp = &(formulation.interpolator());
  interp->interpolate(mesh);

  if(interp->isScalar()){
    InterpolatorScalar& ic = 
      static_cast<InterpolatorScalar&>(*interp); 

    nodalVectorValue = NULL;
    nodalScalarValue = &(ic.getNodeValue()); 
  }
  
  else{
    InterpolatorVector& iv = 
      static_cast<InterpolatorVector&>(*interp); 

    nodalVectorValue = &(iv.getNodeValue());
    nodalScalarValue = NULL; 
 }
  */
}

Solution::~Solution(void){
}

void Solution::write(const string fileName,
		     const string name) const{
  ofstream out;
  out.open(fileName.c_str());

  writeHeader(out);
  writeNodes(out);
  writeElements(out);
  
  writeNodalValues(out, 
		   interp->isScalar(),
		   name);  

  out.close();
}

void Solution::writeHeader(ofstream& out) const{
  out << "$MeshFormat" << endl;
  out << "2.2 0 8" << endl;
  out << "$EndMeshFormat" << endl; 
}

void Solution::writeNodes(ofstream& out) const{
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

void Solution::writeElements(ofstream& out) const{
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
}

void Solution::writeNodalValues(ofstream& out, 
				bool isScalar,
				const string name) const{

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
	out << (*nodalScalarValue)[id] << " ";
      else
	out << (*(*nodalVectorValue)[id])(0) << " "
	    << (*(*nodalVectorValue)[id])(1) << " "
	    << "0 ";
    }

    out << endl;
  }
  out << "$EndElementNodeData" << endl;
}
