#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"
#include "WriterDummy.h"

#include "FormulationPoisson.h"
#include "PoissonSquare.h"
#include "PoissonCircle.h"

using namespace std;

vector<double> fPoisson(Mesh& msh, GroupOfElement& visuDomain, Writer& writer, int order);
vector<double> aPoisson(GroupOfElement& domain, Writer& writer);
fullMatrix<double> l2(fullMatrix<vector<double> >& fem, vector<double>& ana);

int main(int argc, char** argv){
  // Writer //
  WriterMsh writer; 

  // Get Data //
  const unsigned int M        = argc - 3; // Mesh number (without visu)
  const unsigned int maxOrder = atoi(argv[argc - 1]); // Max Order
  Mesh               visu(argv[argc - 2]);
  GroupOfElement     visuDomain = visu.getFromPhysical(9);

  // Analytical Solutions //
  vector<double> ana = aPoisson(visuDomain, writer);

  // Compute FEM //
  fullMatrix<vector<double> > sol(maxOrder, M);

  // Iterate on Meshes
  for(unsigned int i = 0; i < M; i++){
    // Get Mesh
    cout << "** " << argv[1 + i] << endl;
    Mesh msh(argv[1 + i]);

    // Iterate on Orders
    for(unsigned int j = 0; j < maxOrder; j++)
      sol(j, i) = fPoisson(msh, visuDomain, writer, j + 1);
  }
  
  // L2 Error //
  fullMatrix<double> l2Error = l2(sol, ana);
  
  const unsigned int l2Row      = l2Error.size1();
  const unsigned int l2ColMinus = l2Error.size2() - 1;

  cout << "l2 = ..." << endl
       << "    [..." << endl;
  
  for(unsigned int i = 0; i < l2Row; i++){
    cout << "        ";

    for(unsigned int j = 0; j < l2ColMinus; j++)
      cout << scientific << showpos 
	   << l2Error(i, j) << " , ";
    
    cout << scientific << showpos 
	 << l2Error(i, l2ColMinus) << " ; ..." << endl;
  }

  cout << "    ];" << endl;
  
  return 0;
}

vector<double> fPoisson(Mesh& msh, GroupOfElement& visuDomain, Writer& writer, int order){
  // FEM Solution
  GroupOfElement domain = msh.getFromPhysical(9);

  FormulationPoisson poisson(domain, order);
  System sysPoisson(poisson);
  
  sysPoisson.fixDof(msh.getFromPhysical(5), 0);
  sysPoisson.fixDof(msh.getFromPhysical(6), 0);
  sysPoisson.fixDof(msh.getFromPhysical(7), 0);
  sysPoisson.fixDof(msh.getFromPhysical(8), 0);
  
  cout << "Poisson (" << order << "): " << sysPoisson.getSize() << endl;

  sysPoisson.assemble();
  sysPoisson.solve();

  Solution solPoisson(sysPoisson, visuDomain);

  stringstream stream;
  stream << "poisson_Mesh" << domain.getNumber() << "_Order" << order;

  solPoisson.write(stream.str(), writer);

  return solPoisson.getNodalScalarValue();
}

vector<double> aPoisson(GroupOfElement& domain, Writer& writer){
  // Analytical Solution
  cout << "Poisson (Ref)" << endl;

  PoissonCircle poisson(domain);
  poisson.write("poissonRef", writer);

  return poisson.getNodalScalarValue();
}

fullMatrix<double> l2(fullMatrix<vector<double> >& fem, vector<double>& ana){
  // Init //
  const unsigned int nOrder = fem.size1();
  const unsigned int nMesh  = fem.size2();
  const unsigned int nNode  = ana.size();

  fullMatrix<double> res(nOrder, nMesh);

  for(unsigned int i = 0; i < nOrder; i++)
    for(unsigned int j = 0; j < nMesh; j++)
      res(i , j) = 0;

  // Norm of Analytic Solution //
  double anaNorm = 0;
  for(unsigned int k = 0; k < nNode; k++)
    anaNorm += ana[k];// * ana[k];
  
  //anaNorm = sqrt(anaNorm);
  
  // Norm of FEM Error //
  for(unsigned int i = 0; i < nOrder; i++){
    for(unsigned int j = 0; j < nMesh; j++){
      for(unsigned int k = 0; k < nNode; k++)
	res(i, j) += (ana[k] - fem(i, j)[k]);// * (ana[k] - fem(i, j)[k]);
      
      res(i, j) = res(i, j) / anaNorm;//sqrt(res(i, j)) / anaNorm;
    }
  }

  return res;
}
