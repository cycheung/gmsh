#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *loader = 0;

std::string modelName="";
PromptUser *OL = new PromptUser("onelab");
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("gmsh");
InterfacedElast *mySolver = new InterfacedElast("elast");

int main(int argc, char *argv[]){
  bool analyzeOnly=false;
  std::string sockName = "";
  int modelNumber=0;

  getOptions(argc, argv, modelNumber, analyzeOnly, modelName, sockName);

  if (sockName.size())
    loader = new onelab::remoteNetworkClient("loader", sockName);
  Msg::InitializeOnelab("elmer",""); // as a localnetworkclient

  if (loader)
    std::cout << "ONELAB: " << Msg::Synchronize_Down(loader) << " parameters downloaded" << std::endl;

  if (analyzeOnly)
    analyze();
  else
    compute();

  if (loader){
    std::cout << "ONELAB: " << Msg::Synchronize_Up(loader) << " parameters uploaded" << std::endl;
    delete loader;
  }

  Msg::FinalizeOnelab();
}


int analyze(){
  checkIfPresent(modelName+".geo");
  myMesher->analyze("",modelName);
  checkIfPresent(modelName+".dat_onelab");
  mySolver->analyze("",modelName+".dat");
  std::cout << OL->showParamSpace() << std::endl;
  return 1;
}

int compute(){
  newStep();   
  analyze();

  myMesher->run("-2", modelName);
  checkIfModified(modelName+".msh");
  OL->setString("Gmsh/MshFileName", modelName+".msh");

  mySolver->convert(modelName+".dat");
  checkIfModified(modelName+".dat");
  mySolver->run("",modelName);

  std::cout << "Simulation completed successfully" << std::endl;

  GmshDisplay(loader,modelName,OL->getChoices("elast/9OutputFiles"));

  return 1;
}
