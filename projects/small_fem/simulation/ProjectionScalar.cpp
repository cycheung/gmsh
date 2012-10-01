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
		      fullVector<vector<double> >& ana);


unsigned int POW = 1;

double f(fullVector<double>& xyz){
  return pow(xyz(0), POW) + pow(xyz(1), POW);
}

int main(int argc, char** argv){
  // Writer //
  WriterDummy writer;  

  // Get Mesh //
  Mesh msh(argv[1]);
  Mesh visuMesh(argv[2]);

  // Get Domain //
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement visu   = visuMesh.getFromPhysical(7);

  // Get Max Ordre //
  const unsigned int maxOrder = atoi(argv[3]);

  // Max Power //
  const unsigned int P = atoi(argv[4]); // Get Max Power
  POW                  = 1;             // Init POW


  // Iterate On Powers //
  fullVector<vector<double> > real(P);
  fullMatrix<vector<double> > sol(maxOrder, P);

  for(unsigned int p = 0; p < P; p++){

    // Real Function 
    real(p) = ana(visu, f, writer);

    // FEM Solution
    // Iterate on Orders
    for(unsigned int j = 0; j < maxOrder; j++)
      sol(j, p) = fem(domain, visu, f, writer, j + 1);

    // Next power
    POW++;
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

  stream << "projection(" << order << ", " << POW << ")";
  cout   << stream.str()  << ": " << sysProj.getSize() << endl;

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
  stream << "projection(Ref, " << POW << ")";
  cout   << stream.str() << endl;

  Solution projection(f, domain);
  projection.write(stream.str(), writer);

  return projection.getNodalScalarValue();
}


fullMatrix<double> l2(fullMatrix<vector<double> >& fem, 
		      fullVector<vector<double> >& ana){
  // Init //
  const unsigned int nOrder = fem.size1();
  const unsigned int nPower = fem.size2();
  const unsigned int nNode  = ana(0).size();

  fullMatrix<double> res(nOrder, nPower);

  for(unsigned int i = 0; i < nOrder; i++)
    for(unsigned int j = 0; j < nPower; j++)
      res(i , j) = 0;

  // Iterate on Powers //
  for(unsigned int p = 0; p < nPower; p++){

    // Norm of Analytic Solution
    double anaNorm = 0;
    for(unsigned int k = 0; k < nNode; k++)
      anaNorm += ana(p)[k] * ana(p)[k];
  
    anaNorm = sqrt(anaNorm);
  
    // Norm of FEM Error
    for(unsigned int j = 0; j < nOrder; j++){
      for(unsigned int k = 0; k < nNode; k++)
	res(j, p) += (ana(p)[k] - fem(j, p)[k]) * (ana(p)[k] - fem(j, p)[k]);
      
      res(j, p) = sqrt(res(j, p)) / anaNorm;
    }
  }

  return res;
}
