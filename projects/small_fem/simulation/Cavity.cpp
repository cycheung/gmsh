#include <cmath>
#include <cstdio>
#include <iostream>

#include "Mesh.h"
#include "SystemEigen.h"

#include "WriterMsh.h"
#include "Interpolator.h"

#include "FormulationEigenFrequency.h"

#include "SmallFem.h"

using namespace std;

fullVector<double> fDirichlet(fullVector<double>& xyz){
  fullVector<double> f(3);

  f(0) = 0;
  f(1) = 0;
  f(2) = 0;

  return f;
}

void compute(const Options& option){
  cout << option.toString() << endl;
  // Get Domain //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement border = msh.getFromPhysical(5);

  // Writer //
  WriterMsh writer;

  // Get Parameters //
  const size_t order = atoi(option.getValue("-o")[0].c_str());
  const size_t nWave = atoi(option.getValue("-n")[0].c_str());

  // EigenFrequency //
  FormulationEigenFrequency cavity(domain, order);
  SystemEigen sysCavity(cavity);

  //sysCavity.fixCoef(msh.getFromPhysical(5), 0);
  sysCavity.dirichlet(border, fDirichlet);

  cout << "Cavity: " << sysCavity.getSize() << endl;

  cout << "Assembling..." << endl << flush;
  sysCavity.assemble();

  cout << "Solving..." << endl << flush;
  sysCavity.setNumberOfEigenValues(nWave);
  sysCavity.solve();

  // Display //
  const size_t nEigenValue =
    sysCavity.getEigenValuesNumber();

  const vector<complex<double> >& eigenValue =
    sysCavity.getEigenValues();

  cout << "Number of found Eigenvalues: " << nEigenValue
       << endl;

  cout << endl
       << "Number\tEigen Value\tEigen Wave Number" << endl;

  for(size_t i = 0; i < nEigenValue; i++)
    cout << "#" << i + 1        << "\t"
         << eigenValue[i]       << "\t"
         << sqrt(eigenValue[i]) << endl;

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
        sprintf(fileName, "cavity_mode%0*u", nDec, (unsigned int)(i + 1));

        Interpolator intCavity(sysCavity, i, visu);
        intCavity.write(string(fileName), writer);
      }
    }

    else{
      // Without VisuMesh
      for(size_t i = 0; i < nEigenValue; i++){
        sprintf(fileName, "cavity_mode%0*u", nDec, (unsigned int)(i + 1));

        writer.setValues(sysCavity, i);
        writer.write(string(fileName));
      }
    }
  }
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
