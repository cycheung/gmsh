#include <iostream>
#include <sstream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"
#include "WriterVector.h"

#include "FormulationLaplace.h"
#include "FormulationPoisson.h"
#include "FormulationProjection.h"

#include "BasisTest.h"

using namespace std;

int run(int argc, char** argv);
void fLaplace(Mesh& msh, Writer& mWriter);
void fPoisson(Mesh& msh, Mesh& visu, Writer& mWriter, int order);
void fProjection(Mesh& msh, Writer& mWriter);

int main(int argc, char** argv){
  int ret;

  //ret = basisTest(argc, argv);
  ret = run(argc, argv);

  return ret;
}

int run(int argc, char** argv){
  // Writer //
  WriterMsh    mWriter;  
  WriterVector vWriter;
  
  // Get Mesh //
  Mesh msh(argv[1]);
  //cout << msh.toString() << endl;
  
  if(argc == 4){
    Mesh visu(argv[2]);
    int maxOrder = atoi(argv[3]);

    for(int i = 1; i <= maxOrder; i++)
      fPoisson(msh, visu,  mWriter, i);
  }

  //fLaplace(msh, mWriter);
  //fProjection(msh, mWriter);

  return 0;
}

void fLaplace(Mesh& msh, Writer& mWriter){
  GroupOfElement domain = msh.getFromPhysical(7);
  
  FormulationLaplace laplace(domain, 1);
  System sysLaplace(laplace);

  sysLaplace.fixBC(msh.getFromPhysical(6), -1);
  sysLaplace.fixBC(msh.getFromPhysical(5),  2);

  sysLaplace.assemble();
  cout << "Laplace: " << sysLaplace.getSize() << endl;
  sysLaplace.solve();

  Solution solLaplace(sysLaplace);
  solLaplace.write("laplace", mWriter);
}

void fPoisson(Mesh& msh, Mesh& visu, Writer& mWriter, int order){
  GroupOfElement domain = msh.getFromPhysical(9);

  FormulationPoisson poisson(domain, order);
  System sysPoisson(poisson);
  
  sysPoisson.fixBC(msh.getFromPhysical(5), 0);
  sysPoisson.fixBC(msh.getFromPhysical(6), 0);
  sysPoisson.fixBC(msh.getFromPhysical(7), 0);
  sysPoisson.fixBC(msh.getFromPhysical(8), 0);
  
  sysPoisson.assemble();
  cout << "Poisson: " << sysPoisson.getSize() << endl;
  sysPoisson.solve();

  GroupOfElement visuDomain = visu.getFromPhysical(9);
  Solution solPoisson(sysPoisson, visuDomain);

  stringstream stream;
  stream << "poisson" << order;

  solPoisson.write(stream.str(), mWriter);
}

void fProjection(Mesh& msh, Writer& mWriter){
  GroupOfElement domain = msh.getFromPhysical(7);

  fullVector<double> f(3); 
  f(0) =  0.5; 
  f(1) = -1; 
  f(2) =  0; // Vector to project

  FormulationProjection projection(domain, f);
  System sysProj(projection);

  sysProj.assemble();
  cout << "Projection: " << sysProj.getSize() << endl;
  sysProj.solve();

  Solution solProj(sysProj);
  solProj.write("projection", mWriter);
}
