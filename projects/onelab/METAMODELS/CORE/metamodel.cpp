#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *loader = 0;

bool analyzeOnly=false;
std::string modelName = "";

//Onelab clients of the metamodel
PromptUser *OL = new PromptUser("shell"); 
//EncapsulatedTest *myTest = new EncapsulatedTest("../../remote");
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("gmsh");
EncapsulatedGetdp *mySolver = new EncapsulatedGetdp("getdp");

int main(int argc, char *argv[]){
  int modelNumber=0;
  std::string sockName = "";

  OL->setVerbosity(0); //default

  getOptions(argc, argv, modelNumber, analyzeOnly, modelName, sockName);
  modelName="test";

  loader = new onelab::remoteNetworkClient("loader", sockName);
  Msg::InitializeOnelab("metamodel",""); // _onelabClient = new onelab::localClient("metamodel");

  if (sockName.size())
    std::cout << "ONELAB: " << Msg::Synchronize_Down(loader) << " parameters downloaded" << std::endl;
 
  std::cout << OL->showParamSpace();

  if (analyzeOnly)
    analyze();
  else{
    switch (modelNumber){
    case 1:
      OL->setVerbosity(4);
      OL->setNumber("Parameters/Geometry/2Val_Rext", 0.5);
      compute();
      break;
    case 0:
    default:
      compute();
      break;
    }
  }

  if (sockName.size()){
    std::cout << "ONELAB: " << Msg::Synchronize_Up(loader) << " parameters uploaded" << std::endl;
    delete loader;
  }

  Msg::FinalizeOnelab();
}

int analyze(){
  myMesher->analyze("",modelName);
  mySolver->analyze("",modelName);
  std::cout << OL->showParamSpace() << std::endl;
  return 1;
}

int compute(){

  newStep();

  //myTest->run("", modelName)

  checkIfPresent(modelName+".geo");
  myMesher->run("-2", modelName);
  checkIfModified(modelName+".msh");
  OL->setString("Gmsh/MshFileName", modelName+".msh");
   
  OL->setString("GetDP/1ModelName", modelName+".pro");
  checkIfPresent(modelName+".pro");
  mySolver->run("-sol MagSta_phi",modelName);
  checkIfModified(modelName+".res");

  mySolver->run("-pos phi",modelName);
  checkIfModified("phi.pos");
  checkIfModified("b_phi.pos");

  std::cout << "Simulation completed successfully" << std::endl;

  OL->setString("GetDP/9Output files","phi.pos");
  GmshDisplay(loader,modelName,OL->getChoices("GetDP/9Output files"));

  std::cout << OL->showParamSpace() << std::endl;	      

  return 1;
}
