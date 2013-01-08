#include <iostream>

#include "Mesh.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);
    
  Mesh msh(argv[1]);
  cout << msh.toString() << endl;

  GmshFinalize();
  return 0;
}
