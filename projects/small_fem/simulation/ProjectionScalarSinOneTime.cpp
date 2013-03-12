#include <iostream>
#include <sstream>
#include <cmath>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Interpolator.h"
#include "WriterMsh.h"
#include "WriterDummy.h"
#include "Exception.h"

#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FormulationProjectionScalar.h"

#include "Timer.h"
#include "Gmsh.h"

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
  return xyz(0) + xyz(1) + xyz(2);/*
    sin(10 * xyz(0)) +
    sin(10 * xyz(1)) +
    sin(10 * xyz(2));*/
}

int main(int argc, char** argv){
  // Timer //
  Timer timer;

  // Init //
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
  vector<double> real = ana(visu, f, writer);
  cout << "## ... Done !" << endl << flush;


  // Compute FEM //
  fullMatrix<vector<double> > sol(1, 1);
  cout << "## Computing FEM Solution on: " << argv[1]
       << " (Order " << order << ") ..." << endl << flush;

  timer.start();
  sol(0, 0) = fem(domain, visu, f, writer, order);
  timer.stop();

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

  // Done //
  GmshFinalize();

  cout << endl
       << "## Time: "
       << timer.time() << " "
       << timer.unit() << endl;

  return 0;
}

vector<double> fem(GroupOfElement& domain, GroupOfElement& visu,
		   double (*f)(fullVector<double>& xyz),
		   Writer& writer, int order){

  stringstream stream;

  Basis* basis  = BasisGenerator::generate(domain.get(0).getType(),
                                           0, order, "hierarchical");

  FunctionSpaceScalar fSpace(domain, *basis);
  FormulationProjectionScalar projection(f, fSpace);
  System sysProj(projection);

  stream << "projection_Mesh" << domain.getNumber() << "_Order" << order;
  cout   << stream.str()      << ": " << sysProj.getSize() << endl << flush;

  sysProj.assemble();
  sysProj.solve();

  Interpolator intProj(sysProj, visu);

  intProj.write(stream.str(), writer);

  return intProj.getNodalScalarValue();
}

vector<double> ana(GroupOfElement& domain,
		   double (*f)(fullVector<double>& xyz),
		   Writer& writer){

  stringstream stream;

  // Analytical Solution
  stream << "projection_Ref";
  cout   << stream.str() << endl << flush;

  Interpolator projection(f, domain);
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
