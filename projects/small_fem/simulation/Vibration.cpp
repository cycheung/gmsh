#include <cmath>
#include <cstdio>
#include <iostream>

#include "Mesh.h"
#include "SystemEigen.h"

#include "WriterMsh.h"
#include "WriterDummy.h"
#include "Interpolator.h"

#include "FormulationVibration.h"

#include "Gmsh.h"

using namespace std;

double fDir(fullVector<double>& xyz){
  return 0;
}

int main(int argc, char** argv){
  // Init //
  GmshInitialize(argc, argv);
  GmshSetOption("General", "Terminal", 1.);

  // Get Domain //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);

  // Get Some Data //
  WriterMsh writer;

  const unsigned int order = atoi(argv[2]);
  const unsigned int nWave = atoi(argv[3]);

  // Vibration //
  FormulationVibration vibration(domain, order);
  SystemEigen sysVibration(vibration);

  //sysVibration.fixCoef(msh.getFromPhysical(5), 0);

  cout << "Vibration: " << sysVibration.getSize() << endl;

  cout << "Assembling..." << endl << flush;
  sysVibration.assemble();

  cout << "Solving..." << endl << flush;
  sysVibration.setNumberOfEigenValues(nWave);
  sysVibration.solve();

  // Display //
  const unsigned int nEigenValue =
    sysVibration.getEigenValuesNumber();

  const vector<complex<double> >& eigenValue =
    sysVibration.getEigenValues();

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
	      "%s%0*d", "vibration_mode", dec, i + 1);

      Interpolator intVibration(sysVibration, i, visu);
      intVibration.write(string(fileName), writer);
    }
  }

  else{
    // Without VisuMesh
    for(unsigned int i = 0; i < nEigenValue; i++){
      sprintf(fileName,
	      "%s%0*d", "vibration_mode", dec, i + 1);

      writer.setValues(sysVibration, i);
      writer.write(string(fileName));
    }
  }

  // Done //
  GmshFinalize();
  return 0;
}
