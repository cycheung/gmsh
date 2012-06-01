#include "Solution.h"
#include "Exception.h"
#include "Formulation.h"
#include "InterpolatorScalar.h"
#include "InterpolatorVector.h"

using namespace std;

Solution::Solution(const Mesh& mesh, const Formulation& formulation){
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
    nodalScalarValue = ic.getNodeValue(); 
  }
  
  else{
    InterpolatorVector& iv = 
      static_cast<InterpolatorVector&>(*interp); 

    nodalVectorValue = iv.getNodeValue();
    nodalScalarValue = NULL; 
 }
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

  for(int i = 0; i < N; i++)
    out << (*node)[i]->getId() + 1 << " " 
	<< (*node)[i]->getX()      << " "
	<< (*node)[i]->getY()      << " "
	<< (*node)[i]->getZ()      << endl;
  
  out << "$EndNodes" << endl;
}

void Solution::writeElements(ofstream& out) const{
  out << "$Elements" << endl;
  out << E << endl;
  
  for(int i = 0; i < E; i++){
    out << (*element)[i]->getId() + 1 << " " << "2 2 1 1" << " ";
    
    const int M = (*element)[i]->nEntity();
    for(int j = 0; j < M; j++)
      out << (*element)[i]->getEntity(j).getId() + 1 << " ";
    out << endl;
  }

  out << "$EndElements" << endl;
}

void Solution::writeNodalValues(ofstream& out, 
				bool isScalar,
				const string name) const{
  out << "$ElementNodeData"   << endl
      << "1"                  << endl
      << "\"" << name << "\"" << endl
      << "1"                  << endl
      << "0"                  << endl
      << "3"                  << endl
      << "0"                  << endl;
    
  if(isScalar)
      out << "1" << endl;
  else
      out << "3" << endl;
    
    out << E << endl;
  
  for(int i = 0; i < E; i++){
    out << (*element)[i]->getId() + 1 << " " << "3" << " ";
    
    const int M = (*element)[i]->nEntity();
    for(int j = 0; j < M; j++){
      const int id = (*element)[i]->getEntity(j).getId();

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
