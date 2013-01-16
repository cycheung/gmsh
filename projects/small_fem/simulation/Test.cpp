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

  fullMatrix<double> m(3, 3);
  Jacobian j(domain, m);

  j.computeJacobians();
  j.computeInvertJacobians();

  const vector<const pair<const fullMatrix<double>*, double>*>&
    matrix = j.getInvertJacobian(domain.get(17));

  matrix[1]->first->print();
  cout << matrix[1]->second << endl;

  GmshFinalize();
  return 0;
}
