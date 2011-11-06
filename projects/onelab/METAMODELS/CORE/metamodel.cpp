#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

std::vector<double> variables(5);
std::string modelName="test";
PromptUser *promptMe = new PromptUser("shell");
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("gmsh");
EncapsulatedGetdp *mySolver = new EncapsulatedGetdp("getdp");


int metamodel(int modelNumber){
  switch (modelNumber){
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
  
  printf("FHF server has %d Clients\n", onelab::server::instance()->getNumClients());


  myMesher->analyze(modelName); //populate parameterspace
  promptMe->menu("About to mesh",modelName,myMesher->getName());
  myMesher->run("-2", modelName);

  mySolver->analyze(modelName); //populate parameterspace
  promptMe->menu("About to solve",modelName,mySolver->getName());
  mySolver->sol("",modelName);

  promptMe->menu("About to post",modelName,mySolver->getName());

  mySolver->pos("",modelName);

  promptMe->menu("About to leave",modelName,"");

  return 1;
}
