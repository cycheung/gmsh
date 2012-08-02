#include "Solution.h"

using namespace std;

Solution::Solution(const System& system){
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
