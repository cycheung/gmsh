#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

PromptUser *OL = new PromptUser("onelab");
EncapsulatedClient *myMesher = new EncapsulatedClient("gmsh","gmsh",".geo");
InterfacedClient *mySolver = new InterfacedClient("elmer","ElmerSolver",".sif");

void MetaModel::registerClients(){
  registerClient(myMesher);
  registerClient(mySolver);
  myMesher->setFileName(genericNameFromArgs); // setFileName("cryo"); 
  mySolver->setFileName(genericNameFromArgs); 
}

void MetaModel::analyze(){
  Msg::Info("Metamodel::analyze <%s>",getName().c_str());
  myMesher->analyze(); 
  mySolver->analyze(); 
  std::cout << OL->showParamSpace() << std::endl;
}

void MetaModel::compute(){
  Msg::Info("Metamodel::compute <%s>",getName().c_str());
  newStep();   
 
  myMesher->compute();

  mySolver->setLineOptions("");

  mySolver->convert();

  mySolver->compute();

  //GmshDisplay(Msg::loader,genericNameFromArgs,Msg::GetValue("elmer/9OutputFiles"));
  GmshDisplay(Msg::loader,genericNameFromArgs,Msg::GetOnelabChoices("elmer/9OutputFiles"));

  std::cout << "Simulation completed successfully" << std::endl;
}
