#include <iostream>

#include "Mesh.h"
#include "SystemEigen.h"

#include "FormulationEigenFrequencyScalar.h"
#include "FormulationEigenFrequencyVector.h"

#include "SmallFem.h"

using namespace std;

fullVector<double> fDirichletVec(fullVector<double>& xyz){
  fullVector<double> f(3);

  f(0) = 0;
  f(1) = 0;
  f(2) = 0;

  return f;
}

double fDirichletScal(fullVector<double>& xyz){
  return 0;
}

void compute(const Options& option){
  // Get Domain //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement border = msh.getFromPhysical(5);

  // Get Parameters //
  const size_t order = atoi(option.getValue("-o")[0].c_str());
  const size_t nWave = atoi(option.getValue("-n")[0].c_str());

  // Chose write formulation for Eigenvalues and boundary condition //
  Formulation* eig = NULL;
  SystemEigen* sys = NULL;

  if(option.getValue("-type")[0].compare("vector") == 0){
    eig = new FormulationEigenFrequencyVector(domain, order);
    sys = new SystemEigen(*eig);

    sys->dirichlet(border, fDirichletVec);
    cout << "Vectorial ";
  }

  else if(option.getValue("-type")[0].compare("scalar") == 0){
    eig = new FormulationEigenFrequencyScalar(domain, order);
    sys = new SystemEigen(*eig);

    sys->dirichlet(border, fDirichletScal);
    cout << "Scalar ";
  }

  else
    throw Exception("No -type given");

  cout << "Eigenvalues problem: " << sys->getSize() << endl;

  // Assemble and Solve //
  cout << "Assembling..." << endl << flush;
  sys->assemble();

  cout << "Solving..." << endl << flush;
  sys->setNumberOfEigenValues(nWave);
  sys->solve();

  // Display //
  const size_t nEigenValue =
    sys->getEigenValuesNumber();

  const vector<complex<double> >& eigenValue =
    sys->getEigenValues();

  cout << "Number of found Eigenvalues: " << nEigenValue
       << endl;

  cout << endl
       << "Number\tEigen Value" << endl;

  for(size_t i = 0; i < nEigenValue; i++)
    cout << "#" << i + 1  << "\t"
         << eigenValue[i] << endl;

  // Write Sol //
  if(!option.getValue("-nopos").size()){
    FEMSolution feSol;
    sys->addSolution(feSol);
    feSol.write("eigen_mode");
  }

  // Clean //
  delete sys;
  delete eig;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-nopos,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
