#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *loader = 0;

std::string modelName="";
PromptUser *OL = new PromptUser("onelab");
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("gmsh");
InterfacedElmer *mySolver = new InterfacedElmer("elmer");

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
  checkIfPresent(modelName+".sif_onelab");
  mySolver->analyze("",modelName);
  std::cout << OL->showParamSpace() << std::endl;
  return 1;
}

int compute(){
  newStep();   
  analyze();

  myMesher->run("-2", modelName);
  checkIfModified(modelName+".msh");
  OL->setString("Gmsh/MshFileName", modelName+".msh");

  std::string cmd="ElmerGrid 14 2 " + modelName + ".msh -out " + OL->getString("elmer/MshDirName");
  systemCall(cmd); // no server access needed
  checkIfModified(OL->getString("elmer/MshDirName")+"/mesh.header");

  mySolver->convert(modelName);
  checkIfModified(modelName+".sif");
  mySolver->run("",modelName);

  std::cout << "Simulation completed successfully" << std::endl;

  GmshDisplay(loader,modelName,OL->getChoices("elmer/OutputFiles"));

  return 1;
}
