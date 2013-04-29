#include <iostream>

#include "TetReferenceSpace.h"
#include "Gmsh.h"

using namespace std;


int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  TetReferenceSpace ref;
  cout << ref.toLatex() << endl;

  GmshFinalize();
  return 0;
}
