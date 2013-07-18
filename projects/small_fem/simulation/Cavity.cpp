#include <cmath>
#include <cstdio>
#include <iostream>

#include "Mesh.h"
#include "SystemEigen.h"

#include "WriterMsh.h"
#include "Interpolator.h"

#include "FormulationEigenFrequency.h"

#include "Gmsh.h"

using namespace std;

fullVector<double> fDirichlet(fullVector<double>& xyz){
  fullVector<double> f(3);

  f(0) = 0;
  f(1) = 0;
  f(2) = 0;

  return f;
}

int main(int argc, char** argv){
  // Init //
  GmshInitialize(argc, argv);

  // Get Domain //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement border = msh.getFromPhysical(5);

  // Get Some Data //
  WriterMsh writer;

  const size_t order = atoi(argv[2]);
  const size_t nWave = atoi(argv[3]);

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
  // Number of decimals in nEigenValue
  // Used for '0' pading in sprintf
  char fileName[1024];
  const int nDec = floor(log10(nEigenValue)) + 1;

  if(argc == 5){
    // With VisuMesh
    Mesh           visuMesh(argv[4]);
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

  // Done //
  GmshFinalize();
  return 0;
}
