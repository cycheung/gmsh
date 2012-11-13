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

#include "FormulationProjectionVector.h"

#include "Gmsh.h"

using namespace std;

fullVector<double> f(fullVector<double>& xyz);

vector<fullVector<double> > fem(GroupOfElement& domain, GroupOfElement& visu,
				fullVector<double> (*f)(fullVector<double>& xyz),
				Writer& writer, int order);

vector<fullVector<double> > ana(GroupOfElement& domain,
				fullVector<double> (*f)(fullVector<double>& xyz), 
				Writer& writer);

fullMatrix<double> l2(fullMatrix<vector<fullVector<double> > >& fem, 
		      vector<fullVector<double> >& ana);

double modSquare(fullVector<double>& v);
double modDiffSquare(fullVector<double>& v0, fullVector<double>& v1);


fullVector<double> f(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = sin(10 * xyz(0));
  res(1) = sin(10 * xyz(1));
  res(2) = sin(10 * xyz(2));

  return res;
}

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;  

  // Get Data //
  const unsigned int order = atoi(argv[3]);

  // Get Visu //
  cout << "## Reading Visu Mesh ..." << endl << flush;
  Mesh               visuMsh(argv[2]);
  GroupOfElement     visu     = visuMsh.getFromPhysical(7);
  cout << "## ... Done !" << endl << flush;

  // Get Domain //
  cout << "## Reading Mesh ..." << endl << flush;
  Mesh           msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);
  cout << "## ... Done !" << endl << flush;

  // Real Solutions //
  cout << "## Computing Ref Solution ..." << endl << flush;
  vector<fullVector<double> > real = ana(visu, f, writer);
  cout << "## ... Done !" << endl << flush;


  // Compute FEM //
  fullMatrix<vector<fullVector<double> > > sol(1, 1);
  cout << "## Computing FEM Solution on: " << argv[1] 
       << " (Order " << order << ") ..." << endl << flush;

  sol(0, 0) = fem(domain, visu, f, writer, order);
  cout << "## ... Done !" << endl << flush;

  
  // L2 Error //
  cout << "## Computing Error ..." << endl << flush;
  fullMatrix<double> l2Error = l2(sol, real);
  cout << "## ... Done !" << endl << endl << flush;
  
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
  
  GmshFinalize();
  return 0;
}

vector<fullVector<double> > fem(GroupOfElement& domain, GroupOfElement& visu, 
				fullVector<double> (*f)(fullVector<double>& xyz),
				Writer& writer, int order){

  stringstream stream;
  
  FormulationProjectionVector projection(domain, f, order);
  System sysProj(projection);

  stream << "projection_Mesh" << domain.getNumber() << "_Order" << order;
  cout   << stream.str()      << ": " << sysProj.getSize() << endl << flush;

  sysProj.assemble();
  sysProj.solve();

  Solution solProj(sysProj, visu);

  solProj.write(stream.str(), writer);  

  return solProj.getNodalVectorValue();
}

vector<fullVector<double> > ana(GroupOfElement& domain, 
				fullVector<double> (*f)(fullVector<double>& xyz), 
				Writer& writer){
  
  stringstream stream;

  // Analytical Solution
  stream << "projection_Ref";
  cout   << stream.str() << endl << flush;

  Solution projection(f, domain);
  projection.write(stream.str(), writer);

  return projection.getNodalVectorValue();
}

fullMatrix<double> l2(fullMatrix<vector<fullVector<double> > >& fem, 
		      vector<fullVector<double> >& ana){
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
    anaNorm += modSquare(ana[k]);
  
  anaNorm = sqrt(anaNorm);
  
  // Norm of FEM Error //
  for(unsigned int i = 0; i < nOrder; i++){
    for(unsigned int j = 0; j < nMesh; j++){
      for(unsigned int k = 0; k < nNode; k++)
	res(i, j) += modDiffSquare(ana[k], fem(i, j)[k]);
      
      res(i, j) = sqrt(res(i, j)) / anaNorm;
    }
  }

  return res;
}

double modSquare(fullVector<double>& v){
  const int size = v.size();
  double res = 0;

  for(int i = 0; i < size; i++)
    res += v(i) * v(i);

  return res;
}

double modDiffSquare(fullVector<double>& v0, fullVector<double>& v1){
  const int size = v0.size();

  if(size != v1.size())
    throw Exception("Bad Vector Size (modDiffSquare)");

  double res = 0;

  for(int i = 0; i < size; i++)
    res += (v0(i) - v1(i)) * (v0(i) - v1(i));

  return res;
}
