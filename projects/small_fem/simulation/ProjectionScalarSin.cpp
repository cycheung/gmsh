#include <iostream>
#include <sstream>
#include <cmath>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"
#include "WriterDummy.h"
#include "Exception.h"

#include "FormulationProjectionScalar.h"

using namespace std;

double f(fullVector<double>& xyz);

vector<double> fem(GroupOfElement& domain, GroupOfElement& visu,
		   double (*f)(fullVector<double>& xyz),
		   Writer& writer, int order);

vector<double> ana(GroupOfElement& domain,
		   double (*f)(fullVector<double>& xyz), 
		   Writer& writer);

fullMatrix<double> l2(fullMatrix<vector<double> >& fem, 
		      vector<double>& ana);

double f(fullVector<double>& xyz){
  return sin(10 * xyz(0)) + sin(10 * xyz(1));
}

int main(int argc, char** argv){
  // Writer //
  WriterDummy writer;  

  // Get Data //
  const unsigned int M        = argc - 3; // Mesh number (without visu)
  const unsigned int maxOrder = atoi(argv[argc - 1]); // Max Order
  Mesh               visuMsh(argv[argc - 2]);
  GroupOfElement     visu     = visuMsh.getFromPhysical(7);


  // Real Solutions //
  vector<double> real = ana(visu, f, writer);


  // Compute FEM //
  fullMatrix<vector<double> > sol(maxOrder, M);

  // Iterate on Meshes
  for(unsigned int i = 0; i < M; i++){
    // Get Domain
    cout << "** " << argv[1 + i] << endl;
    Mesh           msh(argv[1 + i]);
    GroupOfElement domain = msh.getFromPhysical(7);
    

    // Iterate on Orders
    for(unsigned int j = 0; j < maxOrder; j++)
      sol(j, i) = fem(domain, visu, f, writer, j + 1);
  }


  // L2 Error //
  fullMatrix<double> l2Error = l2(sol, real);
  
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

vector<double> fem(GroupOfElement& domain, GroupOfElement& visu, 
		   double (*f)(fullVector<double>& xyz),
		   Writer& writer, int order){

  stringstream stream;
  
  FormulationProjectionScalar projection(domain, f, order);
  System sysProj(projection);

  stream << "projection_Mesh" << domain.getNumber() << "_Order" << order;
  cout   << stream.str()      << ": " << sysProj.getSize() << endl;

  sysProj.assemble();
  sysProj.solve();

  Solution solProj(sysProj, visu);

  solProj.write(stream.str(), writer);  

  return solProj.getNodalScalarValue();
}

vector<double> ana(GroupOfElement& domain, 
		   double (*f)(fullVector<double>& xyz), 
		   Writer& writer){
  
  stringstream stream;

  // Analytical Solution
  stream << "projection_Ref";
  cout   << stream.str() << endl;

  Solution projection(f, domain);
  projection.write(stream.str(), writer);

  return projection.getNodalScalarValue();
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
    anaNorm += ana[k] * ana[k];
  
  anaNorm = sqrt(anaNorm);
  
  // Norm of FEM Error //
  for(unsigned int i = 0; i < nOrder; i++){
    for(unsigned int j = 0; j < nMesh; j++){
      for(unsigned int k = 0; k < nNode; k++)
	res(i, j) += (ana[k] - fem(i, j)[k]) * (ana[k] - fem(i, j)[k]);
      
      res(i, j) = sqrt(res(i, j)) / anaNorm;
    }
  }

  return res;
}
