#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

std::vector<double> variables(5);
std::string modelName="test";
PromptUser *OL = new PromptUser("shell");
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("gmsh");
EncapsulatedGetdp *mySolver = new EncapsulatedGetdp("getdp");


int metamodel(int modelNumber){
  switch (modelNumber){
  case 1:
    OL->setNumber("1Geometry/2Val_Rext",0.5);
    simulation();
    break;
  case 0:
  default:
    simulation();
    break;
  }
}


int simulation(){
  // std::string sockname("localhost:");
  // const char *port = strstr(sockname.c_str(), ":");
  // int portno = atoi(port + 1); 
  // std::cout << "sockname=<" << sockname << "> portno=" << portno << std::endl;
  // exit(1);
  
  //printf("FHF server has %d Clients\n", onelab::server::instance()->getNumClients());

  if (OL->getInteractivity()) {
    myMesher->analyze(modelName); //populate parameterspace
    OL->menu("About to mesh",modelName,myMesher->getName());
  }
  myMesher->run("-2", modelName);

  mySolver->analyze(modelName); //populate parameterspace
  if (OL->getInteractivity()) {
    OL->menu("About to solve",modelName,mySolver->getName());
  }
  mySolver->sol("",modelName);

  if (OL->getInteractivity()) {
    OL->menu("About to post",modelName,mySolver->getName());
  }

  mySolver->pos("",modelName);

  if (checkWhetherModified("b_phi.pos"))
    std::cout << "Simulation completed successfully" << std::endl;
}
