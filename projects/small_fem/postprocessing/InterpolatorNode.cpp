#include "InterpolatorNode.h"
#include "Jacobian.h"

using namespace std;

InterpolatorNode::InterpolatorNode(void){
  nodeValue   = NULL;
}

InterpolatorNode::~InterpolatorNode(void){  
  if(nodeValue)
    delete nodeValue;
}

void InterpolatorNode::interpolate(const Mesh& mesh){
  int nNode = mesh.getNbNode();
  
  const vector<Node*> node = mesh.getAllNodes();
 
  nodeValue = new vector<double>(nNode);
  
  for(int i = 0; i < nNode; i++)
    (*nodeValue)[i] = node[i]->getValue();

  gotInterpolation = true;
}
