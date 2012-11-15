#include <cmath>
#include <sstream>
#include <iostream>

#include "Mesh.h"
#include "EigenSystem.h"

#include "WriterMsh.h"
#include "Solution.h"

#include "FormulationVibration.h"

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

  // Vibration //  
  FormulationVibration vibration(domain, order);
  EigenSystem sysVibration(vibration);

  sysVibration.fixCoef(msh.getFromPhysical(5), 0);

  cout << "Vibration: " << sysVibration.getSize() << endl;

  sysVibration.assemble();
  sysVibration.solve(nWave);

  // Display //
  const unsigned int nEigenValue = 
    sysVibration.getEigenValueNumber();
  
  const vector<complex<double> >& eigenValue = 
    sysVibration.getEigenValues();

  cout << endl 
       << "Number\tEigen Value\tEigen Wave Number" << endl;
  
  for(unsigned int i = 0; i < nEigenValue; i++)
    cout << "#" << i + 1        << "\t" 
	 << eigenValue[i]       << "\t"
	 << sqrt(eigenValue[i]) << endl;

  // Write Sol //
  for(unsigned int i = 0; i < nEigenValue; i++){
    stringstream stream;
    stream << "vibration_mode" << i + 1;
    
    Solution solVibration(sysVibration, i);
    solVibration.write(stream.str(), writer);
  }


  // Done //
  GmshFinalize();
  return 0;
}
