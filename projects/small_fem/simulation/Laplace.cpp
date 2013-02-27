#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Interpolator.h"
#include "WriterMsh.h"

#include "FormulationLaplace.h"

#include "Timer.h"
#include "Gmsh.h"

using namespace std;

double f1(fullVector<double>& xyz){
  return -1;
}

double f2(fullVector<double>& xyz){
  return 2;
}

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;

  // Get Mesh //
  Mesh msh(argv[1]);

  // Get Order //
  unsigned int order = atoi(argv[2]);

  // Get Domain //
  GroupOfElement domain    = msh.getFromPhysical(7);
  GroupOfElement boundary0 = msh.getFromPhysical(6);
  GroupOfElement boundary1 = msh.getFromPhysical(5);

  // Laplace //
  Timer timer;
  timer.start();

  FormulationLaplace laplace(domain, order);
  System sysLaplace(laplace);

  sysLaplace.dirichlet(boundary0, f1);
  sysLaplace.dirichlet(boundary1, f2);

  cout << "Laplace (" << order << "): "
       << sysLaplace.getSize() << endl;

  sysLaplace.assemble();
  sysLaplace.solve();

  timer.stop();
  cout << "Time: " << timer.time()
       << " "      << timer.unit()
       << endl;

  if(argc == 4){
    // Interpolated View //
    // Visu Mesh
    Mesh visuMesh(argv[3]);
    GroupOfElement visu = visuMesh.getFromPhysical(7);

    Interpolator interp(sysLaplace, visu);
    interp.write("laplace", writer);
  }

  else{
    // Adaptive View //
    writer.setValues(sysLaplace);
    writer.write("laplace");
  }

  GmshFinalize();
  return 0;
}
