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

int main(int argc, char** argv){
  // Init //
  GmshInitialize(argc, argv);

  // Get Domain //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);

  // Get Some Data //
  WriterMsh writer;

  const unsigned int order = atoi(argv[2]);
  const unsigned int nWave = atoi(argv[3]);

  // EigenFrequency //
  FormulationEigenFrequency cavity(domain, order);
  SystemEigen sysCavity(cavity);

  sysCavity.fixCoef(msh.getFromPhysical(5), 0);

  cout << "Cavity: " << sysCavity.getSize() << endl;

  cout << "Assembling..." << endl << flush;
  sysCavity.assemble();

  cout << "Solving..." << endl << flush;
  sysCavity.setNumberOfEigenValues(nWave);
  sysCavity.solve();

  // Display //
  const unsigned int nEigenValue =
    sysCavity.getEigenValuesNumber();

  const vector<complex<double> >& eigenValue =
    sysCavity.getEigenValues();

  cout << "Number of found Eigenvalues: " << nEigenValue
       << endl;

  cout << endl
       << "Number\tEigen Value\tEigen Wave Number" << endl;

  for(unsigned int i = 0; i < nEigenValue; i++)
    cout << "#" << i + 1        << "\t"
	 << eigenValue[i]       << "\t"
	 << sqrt(eigenValue[i]) << endl;

  // Write Sol //
  // Number of decimals in nEigenValue
  // Used for '0' pading in sprintf
  int tmp = nEigenValue;
  int dec = 0;

  while(tmp != 0){
    dec++;
    tmp /= 10;
  }

  char fileName[1024];

  if(argc == 5){
    // With VisuMesh
    Mesh           visuMesh(argv[4]);
    GroupOfElement visu = visuMesh.getFromPhysical(7);

    for(unsigned int i = 0; i < nEigenValue; i++){
      sprintf(fileName,
	      "%s%0*d", "cavity_mode", dec, i + 1);

      Interpolator intCavity(sysCavity, i, visu);
      intCavity.write(string(fileName), writer);
    }
  }

  else{
    // Without VisuMesh
    for(unsigned int i = 0; i < nEigenValue; i++){
      stringstream stream;
      stream << "cavity_mode" << i + 1;

      writer.setValues(sysCavity, i);
      writer.write(stream.str());
    }
  }

  // Done //
  GmshFinalize();
  return 0;
}
