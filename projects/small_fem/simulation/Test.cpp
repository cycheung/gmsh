#include <iostream>

#include "Mesh.h"
#include "GroupOfElement.h"

#include "Jacobian.h"
#include "fullMatrix.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);
    
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);
  
  fullMatrix<double> m(2, 3);
  Jacobian j(domain, m, 1); 
  
  GmshFinalize();
  return 0;
}
