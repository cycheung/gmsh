#include <sstream>
#include <fstream>

#include <set>
#include <deque>

#include "Mesh.h"
#include "Exception.h"

using namespace std;

Mesh::Mesh(const std::string fileName){ 
  // Get Stream //
  in = new ifstream(fileName.c_str());

  // Goto Nodes Definition //
  skipLine(4);
  *in >> nNode;
  node = new vector<Node*>(nNode);
  
  // Create Nodes //
  for(int i = 0; i < nNode; i++){
    int id;
    double x, y, z;
    *in >> id >> x >> y >> z;
    (*node)[i] = new Node(id - 1, x, y, z);
    // NB: I start counting nodes from 0
  }

  // Goto Node Elements Definitions //
  skipLine(3);
  *in >> nNodeElement;
  skipLine(1);

  // Get Lines //
  nodeElement = NULL;
  nLine = 0;
  bool isLine = false;
  do{
    char ln[512];
    in->getline(ln, 512);
    isLine = parse(ln);
  } 
  while(isLine);

  // Get Other NodeElements //
  for(int i = 1; i < nNodeElement; i++){
    char ln[512];
    in->getline(ln, 512);
    parse(ln);
  }

  delete in;

  // Get EdgeElements //
  nEdgeElement = nNodeElement;
  edgeElement  = new vector<Element*>(nEdgeElement);
  getEdges();
  nEdge = edge->size();
}
 
Mesh::~Mesh(void){
  for(int i = 0; i < nNode; i++)
    delete (*node)[i];
  delete node;
  
  for(int i = 0; i < nNodeElement; i++)
    delete (*nodeElement)[i];
  delete nodeElement;

  for(int i = 0; i < nEdge; i++)
    delete (*edge)[i];
  delete edge;

  for(int i = 0; i < nEdgeElement; i++)
    delete (*edgeElement)[i];
  delete edgeElement;
}

string Mesh::toString(void) const{
  stringstream stream;

  stream << "Mesh Informations" << endl;
  stream << "*****************" << endl << endl;

  stream << "Nodes" << endl;
  stream << "-----" << endl << endl;
  for(int i = 0; i < nNode; i++)
    stream << (*node)[i]->toString() << endl;
  stream << endl;

  stream << "Node Elements" << endl;
  stream << "-------------" << endl << endl;  
  for(int i = 0; i < nNodeElement; i++)
    stream << (*nodeElement)[i]->toString() << endl;
  stream << endl;

  stream << "Edges" << endl;
  stream << "-----" << endl << endl;
  for(int i = 0; i < nEdge; i++)
    stream << (*edge)[i]->toString() << endl;
  stream << endl;

  stream << "Edge Elements" << endl;
  stream << "-------------" << endl << endl;  
  for(int i = 0; i < nEdgeElement; i++)
    stream << (*edgeElement)[i]->toString() << endl;

  return stream.str();
}

void Mesh::skipLine(const int N){
  char ln[512];
  for(int i = 0; i < N; i++)
    in->getline(ln, 512);
}

bool Mesh::parse(const string str){
  istringstream stream(str);
  int id;
  int type;

  stream >> id >> type;

  switch(type){
  case 1: getLine(id, stream); return true;
  case 2: getTri(id, stream);  return false;
  case 3: throw Exception("I don't do Quads !"); return false;
  
  default: return false;
  }
}

void Mesh::getTri(int id, istringstream& stream){
  // Is Element allocated ? //
  if(!nodeElement){
    nNodeElement -= nLine;
    nodeElement = new vector<Element*>(nNodeElement);
  }

  // New Element
  id -= 1 + nLine;
  (*nodeElement)[id] = new Element(id, 30);

  // Get Node
  int nodeId[3];
  stream >> nodeId[0] >> nodeId[0] >> nodeId[0] >> nodeId[0];
  stream >> nodeId[1];
  stream >> nodeId[2];
 
  (*nodeElement)[id]->node->at(0) = (*node)[nodeId[0] - 1];
  (*nodeElement)[id]->node->at(1) = (*node)[nodeId[1] - 1];
  (*nodeElement)[id]->node->at(2) = (*node)[nodeId[2] - 1];

  (*nodeElement)[id]->entity->at(0) = (*node)[nodeId[0] - 1];
  (*nodeElement)[id]->entity->at(1) = (*node)[nodeId[1] - 1];
  (*nodeElement)[id]->entity->at(2) = (*node)[nodeId[2] - 1];

  (*nodeElement)[id]->orientation->at(0) = 1;
  (*nodeElement)[id]->orientation->at(1) = 1;
  (*nodeElement)[id]->orientation->at(2) = 1;

  (*nodeElement)[id]->buildJacobian();
}

void Mesh::getLine(int id, istringstream& stream){
  int phys;
  int nodeId;
  
  stream >> phys >> phys;
  for(int i = 0; i < 3; i++){
    stream >> nodeId;
    nodeId -= 1;
    if(!(*node)[nodeId]->gotPhysical())
      (*node)[nodeId]->setPhysical(phys);
  }

  nLine += 1;
}

void Mesh::getEdges(void){
  deque<Edge*> stack;
  set<edgeTriplet, EdgeTripletComparator> lookup;
  int edgeId = 0;

  for(int i = 0; i < nNodeElement; i++){
    (*edgeElement)[i] = new Element(i, 31); // New Edge Element

    // Get Geometric informations
    const vector<Node*>& nodes = (*nodeElement)[i]->getAllNodes();
    const int N = nodes.size();

    for(int j = 0; j < N; j++)
      (*edgeElement)[i]->node->at(j) = nodes[j];

    (*edgeElement)[i]->buildJacobian();

    // Edge Element are build from node Element
    const vector<Entity*>& n = (*nodeElement)[i]->getAllEntities();
    const int E = n.size();

    // Get *Unique* Edges //
    for(int j = 0; j < E; j++){
      int k = (j + 1) % 3; 
      struct edgeTriplet idUnswap;
      struct edgeTriplet idSwap;
      
      idUnswap.origin = n[j]->getId(); idUnswap.end = n[k]->getId();
      idSwap.origin   = n[k]->getId(); idSwap.end   = n[j]->getId();
      
      pair<set<edgeTriplet>::iterator, bool>
	notInLookUpUnswap = lookup.insert(idUnswap);
      pair<set<edgeTriplet>::iterator, bool>
	notInLookUpSwap   = lookup.insert(idSwap);

      // If New Edge
      if(notInLookUpUnswap.second && 
	 notInLookUpSwap.second){
	struct edgeTriplet* t0;
	struct edgeTriplet* t1;
	t0 = const_cast<struct edgeTriplet*>(&(*notInLookUpUnswap.first));
	t1 = const_cast<struct edgeTriplet*>(&(*notInLookUpSwap.first));
	
	// NOTE
	// const_casts are OK:
	// We change only the 'edge' of our edgeTriplet
	// and 'egde' doesn't change the ordering of the set
	 
	t0->edge = new Edge(edgeId, *((Node*)n[j]), *((Node*)n[k]));
	t1->edge = t0->edge;
	edgeId++;
	
	(*edgeElement)[i]->entity->at(j) = t0->edge;
	(*edgeElement)[i]->orientation->at(j) = 1;

	stack.push_back(t0->edge);
      }
      
      // If Edge Found
      if(!notInLookUpUnswap.second){
	(*edgeElement)[i]->entity->at(j) = notInLookUpUnswap.first->edge;	
	(*edgeElement)[i]->orientation->at(j) = 1;
      }

      if(!notInLookUpSwap.second){
	(*edgeElement)[i]->entity->at(j) = notInLookUpSwap.first->edge;
  	(*edgeElement)[i]->orientation->at(j) = -1; //Edge goes wrong way
      } 
    }
  }

  // Store Edges //
  deque<Edge*>::iterator ie   = stack.begin();
  deque<Edge*>::iterator eend = stack.end();
  
  edge = new vector<Edge*>(ie, eend);
}
