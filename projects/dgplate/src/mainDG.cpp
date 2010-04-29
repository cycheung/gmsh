#include "Gmsh.h"
#include "GmshMessage.h"
#include "Context.h"
#include "GModel.h"
#include "LuaBindings.h"

#include "DgC0PlateSolver.h"
#include "PView.h"
#include "PViewData.h"
#include "groupOfElements.h"
#include <iterator>
#include <iostream>
#include "IPState.h"
int main (int argc, char* argv[])
{

  if (argc != 2){
    printf("usage : elasticity input_file_name\n");
    return -1;
  }
  GmshInitialize(argc, argv);
  GmshSetOption("General","Terminal",1.);

  //binding *b = binding::instance();
  //DgC0PlateSolver::registerBindings(b);
  // globals are still present in Gmsh

  // instanciate a solver
  DgC0PlateSolver mySolver (1000);

  // read some input file
  mySolver.readInputFile(argv[1]);

  // create the interface element after reading the mesh file put this somewhere else ??
  mySolver.createInterfaceElement();
  // There is a bug in the function see modifGMSH
  //mySolver.createInterfaceElement_2();

  // solve the problem
  mySolver.solve();

  PView *pv = mySolver.buildDisplacementView("displacement");
  if(pv!=NULL) pv->getData()->writeMSH("disp.msh", false);
  delete pv;
  mySolver.buildVonMisesView("stressVM.msh");
  //pv = mySolver.buildElasticEnergyView("elastic energy"); //error ??
  //pv->getData()->writeMSH("energ.msh", false);
  //delete pv;

  // stop gmsh
  GmshFinalize();

}
