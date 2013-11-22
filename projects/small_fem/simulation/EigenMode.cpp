#include <cstdio>
#include <iostream>

#include "Mesh.h"
#include "SystemEigen.h"

#include "WriterMsh.h"
#include "Interpolator.h"

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

  // Writer //
  WriterMsh writer;

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
    // Number of decimals in nEigenValue
    // Used for '0' pading in sprintf
    char fileName[1024];
    const int nDec = floor(log10(nEigenValue)) + 1;

    if(option.getValue("-interp").size()){
      // With VisuMesh
      Mesh           visuMesh(option.getValue("-interp")[0]);
      GroupOfElement visu = visuMesh.getFromPhysical(7);

      for(size_t i = 0; i < nEigenValue; i++){
        sprintf(fileName, "eigen_mode%0*u", nDec, (unsigned int)(i + 1));

        Interpolator intCavity(*sys, i, visu);
        intCavity.write(string(fileName), writer);
      }
    }

    else{
      // Without VisuMesh
      for(size_t i = 0; i < nEigenValue; i++){
        sprintf(fileName, "eigen_mode%0*u", nDec, (unsigned int)(i + 1));

        writer.setValues(*sys, i);
        writer.write(string(fileName));
      }
    }
  }

  // Clean //
  delete sys;
  delete eig;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-n,-nopos,-interp,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
