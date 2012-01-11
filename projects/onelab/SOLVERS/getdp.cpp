#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *loader = 0;

std::string modelName = "";
PromptUser *OL = new PromptUser("shell"); 
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("mygmsh");
EncapsulatedGetdp *mySolver = new EncapsulatedGetdp("getdp");

int main(int argc, char *argv[]){
  bool analyzeOnly=false;
  std::string sockName = "";
  int modelNumber=0;

  getOptions(argc, argv, modelNumber, analyzeOnly, modelName, sockName);

  Msg::InitializeOnelab("onelab",""); // _onelabClient = new onelab::localClient("metamodel");

  std::cout << OL->showParamSpace() << std::endl;
  if(sockName.size())
    loader = new onelab::remoteNetworkClient("loader", sockName);

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
  std::cout << OL->showParamSpace() << std::endl;
  Msg::FinalizeOnelab();
}

int analyze(){
  checkIfPresent(modelName+".geo");
  myMesher->analyze("",modelName);
  checkIfPresent(modelName+".pro");
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
  OL->setString("GetDP/1ModelName", modelName+".pro");

  mySolver->run("-sol " + OL->getString("GetDP/1ResolutionChoices"),modelName);
  checkIfModified(modelName+".res");

  mySolver->run("-pos " + OL->getString("GetDP/2PostOperationChoices"),modelName);

  std::cout << "Simulation completed successfully" << std::endl;

  if(loader){
    std::vector<std::string> choices = OL->getChoices("GetDP/9Output files");
    for(unsigned int i = 0; i < choices.size(); i++){
      std::string filename = choices[i];
      checkIfModified(choices[i]);
      loader->sendMergeFileRequest(choices[i]);
    }
  }
  return 1;
}
