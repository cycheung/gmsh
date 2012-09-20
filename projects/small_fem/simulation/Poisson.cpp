#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationPoisson.h"
#include "PoissonSquare.h"

using namespace std;

vector<double> fPoisson(Mesh& msh, Mesh& visu, Writer& writer, int order);
vector<double> aPoisson(Mesh& msh, Writer& writer);
vector<double> l2(vector<vector<double> >& v);

int main(int argc, char** argv){
  // Writer //
  WriterMsh writer; 
  
  // Get Mesh //
  Mesh msh(argv[1]);
  Mesh visu(argv[2]);

  // Compute FEM //
  int maxOrder = atoi(argv[3]);
  vector<vector<double> > sol(maxOrder + 1);

  for(int i = 1; i <= maxOrder; i++)
    sol[i - 1] = fPoisson(msh, visu,  writer, i);

  // Compute Analytical //
  sol[maxOrder] = aPoisson(visu, writer);

  // L2 Error //
  vector<double> l2Error = l2(sol);

  for(unsigned int i = 0; i < l2Error.size(); i++)
    cout << i + 1 << ": " << l2Error[i] << endl;

  return 0;
}

vector<double> fPoisson(Mesh& msh, Mesh& visu, Writer& writer, int order){
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

  solPoisson.write(stream.str(), writer);

  return solPoisson.getNodalScalarValue();
}

vector<double> aPoisson(Mesh& msh, Writer& writer){
  // Analytical Solution
  GroupOfElement domain = msh.getFromPhysical(9);

  cout << "Poisson (Ref)" << endl;
  PoissonSquare poisson(domain);
  poisson.write("poissonRef", writer);

  return poisson.getNodalScalarValue();
}

vector<double> l2(vector<vector<double> >& v){
  unsigned int size      = v.size();
  unsigned int sizeMinus = size - 1;
  unsigned int node      = v[0].size();

  vector<double> res(sizeMinus, 0);

  for(unsigned int i = 0; i < sizeMinus; i++){
    for(unsigned int j = 0; j < node; j++)
      res[i] += (v[i][j] - v[sizeMinus][j]) * (v[i][j] - v[sizeMinus][j]);
    
    res[i] = sqrt(res[i]);
  }

  return res;
}
