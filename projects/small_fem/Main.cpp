#include <iostream>
#include <string>
#include <vector>
#include "Mesh.h"
#include "fullMatrix.h"
#include "FormulationLaplace.h"
#include "FormulationProjection.h"
#include "System.h"
#include "Solution.h"

#include "GModel.h"
#include "MElement.h"

using namespace std;

int main(int argc, char** argv){
  GModel msh("SmallFEM");
  msh.readMSH(string(argv[1]));

  int N = msh.getNumMeshElements();

  cout << "Test: " << N << endl;

  for(int i = 1; i <= N; i++)
    cout << i << ": " << msh.getMeshElementByTag(i)->getType() << endl;

  /*
  // Get Mesh //
  Mesh msh(argv[1]);
  

  // Laplace //
  FormulationLaplace laplace;
  System sysLaplace(msh.getAllNodeElements(), laplace);

  sysLaplace.assemble();

  sysLaplace.fixBC(5, -2);
  sysLaplace.fixBC(6,  1);

  sysLaplace.solve();

  Solution solLaplace(msh, laplace);
  solLaplace.write("laplace.pos", "laplace");

    
  // Projection //
  fullVector<double> f(2); 
  f(0) = -1; f(1) = 1; // Vector to project
  
  FormulationProjection projection(f);
  System sysProj(msh.getAllEdgeElements(), projection);

  sysProj.solve();
  
  Solution solProj(msh, projection);
  solProj.write("projection.pos", "projection");
  */      
  
  return 0;
}

