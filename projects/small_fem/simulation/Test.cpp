#include <iostream>

#include "Mesh.h"
#include "WriterMsh.h"

#include "TetReferenceSpace.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);
    
  Mesh msh(argv[1]);
  //cout << msh.toString() << endl;

  GroupOfElement goe = msh.getFromPhysical(7);

  for(unsigned int i = 0; i < goe.getNumber(); i++){
    const MElement& element = goe.get(i);
    MElement& ce = const_cast<MElement&>(element);

    cout << "Element #" << element.getNum() << endl;
    
    for(int e = 0; e < ce.getNumEdges(); e++)
      cout << " -- Edge: " 
	   << ce.getEdge(e).getVertex(0)->getNum() - 1
	   << ", "
	   << ce.getEdge(e).getVertex(1)->getNum() - 1
	   << endl;
    
    for(int f = 0; f < ce.getNumFaces(); f++)
      cout << " -- Face: " 
	   << ce.getFace(f).getVertex(0)->getNum() - 1
	   << ", "
	   << ce.getFace(f).getVertex(1)->getNum() - 1
	   << ", "
	   << ce.getFace(f).getVertex(2)->getNum() - 1
	   << endl;    
  }
    
  //TetReferenceSpace ref;
  //cout << ref.toString() << endl;

  GmshFinalize();
  return 0;
}
