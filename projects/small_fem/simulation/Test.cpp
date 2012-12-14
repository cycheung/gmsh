#include <iostream>

#include "Mesh.h"
#include "WriterMsh.h"

#include "LineReferenceSpace.h"

#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);
  
  LineReferenceSpace a;
  cout << a.toLatex() << endl;
  cout << a.toString() << endl;  

  GmshFinalize();
  return 0;
}
