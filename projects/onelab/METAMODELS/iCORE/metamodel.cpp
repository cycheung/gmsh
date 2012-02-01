#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

std::vector<double> variables(5);
std::string modelName="core";
PromptUser *OL = new PromptUser("onelab");
InterfacedGmsh *mesher = new InterfacedGmsh("gmsh");
InterfacedGetdp *solver = new InterfacedGetdp("getdp");


int metamodel(int modelNumber){
  switch (modelNumber){
  case 1:
    OL->setNumber("FormulationA",0); 
    simulation();
    break;
  case 0:
  default:
    simulation();
    break;
  }
}

int simulation(){ 
  mesher->analyze(modelName);
  OL->menu("After meshing",modelName,mesher->getName());
  mesher->run("-2", modelName);

  solver->analyze(modelName); // populate parameterspace
  OL->menu("About to solve with GETDP",modelName,solver->getName());
  solver->run("-sol MagSta -pos MagSta", modelName);

  OL->menu("After solving with GETDP",modelName,solver->getName());
}
