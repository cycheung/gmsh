#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"
#include "WriterVector.h"

#include "FormulationLaplace.h"
#include "FormulationPoisson.h"
#include "FormulationProjection.h"

#include "PoissonSquare.h"

#include "BasisTest.h"

using namespace std;

int run(int argc, char** argv);
void fLaplace(Mesh& msh, Writer& mWriter);
vector<double> fPoisson(Mesh& msh, Mesh& visu, Writer& mWriter, int order);
vector<double> aPoisson(Mesh& msh, Writer& mWriter);
void fProjection(Mesh& msh, Writer& mWriter);

vector<double> l2(vector<vector<double> >& v);

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

    vector<vector<double> > sol(maxOrder + 1);

    for(int i = 1; i <= maxOrder; i++)
      sol[i - 1] = fPoisson(msh, visu,  vWriter, i);

    sol[maxOrder] = aPoisson(visu, vWriter);

    vector<double> l2Error = l2(sol);

    for(unsigned int i = 0; i < l2Error.size(); i++)
      cout << i + 1 << ": " << l2Error[i] << endl;
  }

  //fLaplace(msh, mWriter);
  //fProjection(msh, mWriter);

  return 0;
}

void fLaplace(Mesh& msh, Writer& mWriter){
  GroupOfElement domain = msh.getFromPhysical(7);
  
  FormulationLaplace laplace(domain, 1);
  System sysLaplace(laplace);

  sysLaplace.fixDof(msh.getFromPhysical(6), -1);
  sysLaplace.fixDof(msh.getFromPhysical(5),  2);

  sysLaplace.assemble();
  cout << "Laplace: " << sysLaplace.getSize() << endl;
  sysLaplace.solve();

  Solution solLaplace(sysLaplace);
  solLaplace.write("laplace", mWriter);
}

vector<double> fPoisson(Mesh& msh, Mesh& visu, Writer& mWriter, int order){
  // FEM Solution
  GroupOfElement domain = msh.getFromPhysical(9);

  FormulationPoisson poisson(domain, order);
  System sysPoisson(poisson);
  
  sysPoisson.fixDof(msh.getFromPhysical(5), 0);
  sysPoisson.fixDof(msh.getFromPhysical(6), 0);
  sysPoisson.fixDof(msh.getFromPhysical(7), 0);
  sysPoisson.fixDof(msh.getFromPhysical(8), 0);
  
  sysPoisson.assemble();
  cout << "Poisson (" << order << "): " << sysPoisson.getSize() << endl;
  //cout << "Function Space:" << endl << poisson.fs().toString() << endl;
  sysPoisson.solve();

  GroupOfElement visuDomain = visu.getFromPhysical(9);
  Solution solPoisson(sysPoisson, visuDomain);

  stringstream stream;
  stream << "poisson" << order;

  solPoisson.write(stream.str(), mWriter);

  return solPoisson.getNodalScalarValue();
}

vector<double> aPoisson(Mesh& msh, Writer& mWriter){
  // Analytical Solution
  GroupOfElement domain = msh.getFromPhysical(9);

  cout << "Poisson (Ref)" << endl;
  PoissonSquare poisson(domain);
  poisson.write("poissonRef", mWriter);

  return poisson.getNodalScalarValue();
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

vector<double> l2(vector<vector<double> >& v){
  unsigned int size      = v.size();
  unsigned int sizeMinus = size - 1;
  unsigned int node      = v[0].size();

  vector<double> res(sizeMinus, 0);

  for(unsigned int i = 0; i < size; i++){
    for(unsigned int j = 0; j < node; j++)
      res[i] += (v[i][j] - v[sizeMinus][j]) * (v[i][j] - v[sizeMinus][j]);
    
    res[i] = sqrt(res[i]);
  }

  return res;
}
