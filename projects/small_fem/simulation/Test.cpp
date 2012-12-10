#include <iostream>

#include "Mesh.h"
#include "WriterMsh.h"

#include "TriReferenceSpace.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);
  
  // Writer //
  WriterMsh writer;
  
  // Get Mesh //
  Mesh msh(argv[1]);
 
  // Get Domain //
  GroupOfElement domain = 
    msh.getFromPhysical(7);

  TriReferenceSpace ref;
  
  cout << ref.toString() << endl << endl;

  unsigned int id = ref.getReferenceSpace(domain.get(1));
  cout << id << endl;
  

  vector<const vector<const vector<unsigned int>*>*> edge = 
    ref.getAllEdge();

  for(unsigned int i = 0; i < edge[id]->size(); i++)
    cout << "Edge " << i << ": (" 
	 << edge[id]->at(i)->at(0) << ", "
	 << edge[id]->at(i)->at(1) << ")" << endl;
  
  GmshFinalize();
  return 0;
}
