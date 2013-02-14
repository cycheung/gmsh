#include <iostream>

#include "Dof.h"
#include "Mesh.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Mesh msh(argv[1]);
  // cout << msh.toString() << endl;

  Dof  d(1, 1);
  Dof* ptr;

  cout << sizeof(d)
       << " | "
       << sizeof(ptr)
       << endl;

  GmshFinalize();
  return 0;
}
