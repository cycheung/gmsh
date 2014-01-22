#include <iostream>
#include <complex>

#include "Mesh.h"
#include "SystemEigen.h"
#include "SystemHelper.h"

#include "FormulationEigenFrequencyScalar.h"
#include "FormulationEigenFrequencyVector.h"

#include "SmallFem.h"

using namespace std;

fullVector<complex<double> > fDirichletVec(fullVector<double>& xyz){
  fullVector<complex<double> > f(3);

  f(0) = complex<double>(0, 0);
  f(1) = complex<double>(0, 0);
  f(2) = complex<double>(0, 0);

  return f;
}

complex<double> fDirichletScal(fullVector<double>& xyz){
  return complex<double>(0, 0);
}

void compute(const Options& option){
  // Get Domain //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement border = msh.getFromPhysical(5);

  // Get Parameters //
  const size_t order = atoi(option.getValue("-o")[1].c_str());
  const size_t nWave = atoi(option.getValue("-n")[1].c_str());

  // Chose write formulation for Eigenvalues and boundary condition //
  Formulation<complex<double> >* eig = NULL;
  SystemEigen* sys = NULL;

  if(option.getValue("-type")[1].compare("vector") == 0){
    eig = new FormulationEigenFrequencyVector(domain, order);
    sys = new SystemEigen(*eig);

    SystemHelper<complex<double> >::dirichlet(*sys, border, fDirichletVec);
    cout << "Vectorial ";
  }

  else if(option.getValue("-type")[1].compare("scalar") == 0){
    eig = new FormulationEigenFrequencyScalar(domain, order);
    sys = new SystemEigen(*eig);

    SystemHelper<complex<double> >::dirichlet(*sys, border, fDirichletScal);
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
  fullVector<complex<double> > eigenValue;
  const size_t nEigenValue = sys->getNComputedSolution();
  sys->getEigenValues(eigenValue);

  cout << "Number of found Eigenvalues: " << nEigenValue
       << endl;

  cout << endl
       << "Number\tEigen Value" << endl;

  for(size_t i = 0; i < nEigenValue; i++)
    cout << "#" << i + 1  << "\t"
         << eigenValue(i) << endl;

  // Write Sol //
  try{
    option.getValue("-nopos");
  }
  catch(...){
    FEMSolution<complex<double> > feSol;
    sys->getSolution(feSol);
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
