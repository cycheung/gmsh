#include <cmath>
#include <cstdio>
#include <iostream>

#include "Mesh.h"
#include "SystemEigen.h"

#include "WriterMsh.h"
#include "Interpolator.h"

#include "FormulationVibration.h"

#include "SmallFem.h"

using namespace std;

double fDir(fullVector<double>& xyz){
  return 0;
}

void compute(const Options& option){
  // Get Domain //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement border = msh.getFromPhysical(5);

  // Writer //
  WriterMsh writer;

  // Get Some Data //
  const size_t order = atoi(option.getValue("-o")[0].c_str());
  const size_t nWave = atoi(option.getValue("-n")[0].c_str());

  // Vibration //
  FormulationVibration vibration(domain, order);
  SystemEigen sysVibration(vibration);

  //sysVibration.fixCoef(msh.getFromPhysical(5), 0);
  sysVibration.dirichlet(border, fDir);

  cout << "Vibration: " << sysVibration.getSize() << endl;

  cout << "Assembling..." << endl << flush;
  sysVibration.assemble();

  cout << "Solving..." << endl << flush;
  sysVibration.setNumberOfEigenValues(nWave);
  sysVibration.solve();

  // Display //
  const size_t nEigenValue =
    sysVibration.getEigenValuesNumber();

  const vector<complex<double> >& eigenValue =
    sysVibration.getEigenValues();

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
        sprintf(fileName, "vibration_mode%0*u", nDec, (unsigned int)(i + 1));

        Interpolator intVibration(sysVibration, i, visu);
        intVibration.write(string(fileName), writer);
      }
    }

    else{
      // Without VisuMesh
      for(size_t i = 0; i < nEigenValue; i++){
        sprintf(fileName, "vibration_mode%0*u", nDec, (unsigned int)(i + 1));

        writer.setValues(sysVibration, i);
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
