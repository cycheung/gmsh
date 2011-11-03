#include <stdlib.h>
#include <string>
#include "OnelabClients.h"


std::string modelName="cryo";
PromptUser *OL = new PromptUser("onelab");
InterfacedElmer *elmer = new InterfacedElmer("elmer");
std::vector<double> variables(5);


int metamodel(int modelNumber){
  switch (modelNumber){
  case 1:
    variables[0] = 77;
    std::cout << "With Tcold =" << variables[0] << " K " << std::endl;
    simulation();
    std::cout << "\nThe optimum time is " << variables[1] << " sec." << std::endl;
    break;
  case 0:
  default:
    simulation();
    break;
  }
}

int simulation(){ 

  elmer->analyze(modelName); // populate parameterspace

  if (OL->getInteractivity())
    OL->menu("About to convert",modelName,elmer->getName());

  elmer->convert(modelName); 

  OL->setNumber("Tcold",variables[0]); // overrides default value
  OL->setNumber("NumStep",50); // overrides default value
  OL->setNumber("TimeStep",0.05); // overrides default value

  if (OL->getInteractivity())
    OL->menu("About to solve with ELMER",modelName,elmer->getName());

  elmer->run("",modelName);

  array data = read_array("tempevol.txt",' '); // temperature over time at probe point
  variables[1] = find_in_array((int)(1./OL->getNumber("TimeStep")),2,data); // temperature after 2s
  
  system("gmsh - solution.pos script.pos.opt &> gmsh.script.log"); // no server access needed

  std::vector<double> v1=extract_column(8,read_array("f1.txt",' '));
  std::vector<double> v2=extract_column(8,read_array("f2.txt",' '));
  std::vector<double> time=extract_column(2,read_array("f1.txt",' '));
  double min=1e200, temp;
  int imin=0;
  
  for (int i=0; i<v1.size(); i++){
    temp = v1[i]+v2[i];
    if (temp<min){ imin=i; min=temp;}
  }
  variables[1] = time[imin];
}
