#include "InterpolatorEdge.h"
#include "Jacobian.h"

using namespace std;

InterpolatorEdge::InterpolatorEdge(const BasisVector& basis){
  this->basis    = basis.getBasis();
  this->bSize    = basis.getSize();
  nNode          = 0;
  nodeValue      = NULL;
  isInterpolated = NULL;
}

InterpolatorEdge::~InterpolatorEdge(void){
  for(int i = 0; i < nNode; i++)
    delete (*nodeValue)[i];

  if(nodeValue)
    delete nodeValue;
  
  if(isInterpolated)
    delete isInterpolated;
}

void InterpolatorEdge::interpolate(const Mesh& mesh){
  msh   = &mesh;
  nNode =  mesh.getNbNode();
  
  nodeValue      = new vector<fullVector<double>*>(nNode);
  isInterpolated = new vector<bool>(nNode);

  for(int i = 0; i < nNode; i++){
    (*nodeValue)[i] = new fullVector<double>(2);

    (*nodeValue)[i]->set(0, 0.0); 
    (*nodeValue)[i]->set(1, 0.0);
    
    (*isInterpolated)[i] = false;
  }
  
  interpolateEdgeElement();
}

void InterpolatorEdge::interpolateEdgeElement(void){
  const vector<Element*>& element = msh->getAllEdgeElements();
  const int E = element.size();

  for(int i = 0; i < E; i++){
    const vector<Node*>&   node   = element[i]->getAllNodes();
    const vector<Entity*>& entity = element[i]->getAllEntities();
    const vector<int>&     orient = element[i]->getAllOrientations();

    const int N = node.size();
    const Jacobian& jac = element[i]->getJacobian();

    for(int j = 0; j < N; j++){
      const int id   = node[j]->getId();
      
      if(!(*isInterpolated)[id]){
	const double x = node[j]->getX();
	const double y = node[j]->getY();
	
	fullVector<double>* vn = (*nodeValue)[id];
	fullVector<double>  uv = jac.invMap(x, y);
	
	for(int k = 0; k < bSize; k++){
	  fullVector<double> vk = jac.grad(basis[k].at(uv(0), uv(1), 0));

	  vk.scale(entity[k]->getValue() * orient[k]);
	  vn->axpy(vk, 1);
	}
	
	(*isInterpolated)[id] = true;
      }
    }
  }

  gotInterpolation = true;
}
