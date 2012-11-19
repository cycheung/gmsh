#include <cmath>
#include <iostream>

#include "Mesh.h"
#include "EigenSystem.h"

#include "WriterMsh.h"
#include "Solution.h"

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
  EigenSystem sysCavity(cavity);

  sysCavity.fixCoef(msh.getFromPhysical(5), 0);

  cout << "Cavity: " << sysCavity.getSize() << endl;

  cout << "Assembling..." << endl << flush;
  sysCavity.assemble();

  cout << "Solving..." << endl << flush;
  sysCavity.solve(nWave);

  // Display //
  const unsigned int nEigenValue = 
    sysCavity.getEigenValueNumber();
  
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
  if(argc == 5){
    // With VisuMesh
    Mesh           visuMesh(argv[4]);
    GroupOfElement visu = visuMesh.getFromPhysical(7);
  
    for(unsigned int i = 0; i < nEigenValue; i++){
      stringstream stream;
      stream << "cavity_mode" << i + 1;
      
      Solution solCavity(sysCavity, i, visu);
      solCavity.write(stream.str(), writer);
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
