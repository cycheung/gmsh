#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

PromptUser *OL = new PromptUser("onelab");
void MetaModel::analyze(){
  Msg::Info("Metamodel::analyze <%s>",getName().c_str());
  simpleCheck();
  //std::cout << OL->showParamSpace() << std::endl;
  std::cout << OL->showClientStatus() << std::endl;
}

void MetaModel::compute(){
  Msg::Info("Metamodel::compute <%s>",getName().c_str());
  simpleCompute();
  GmshDisplay(Msg::loader,Msg::GetOnelabString("Arguments/FileName"),Msg::GetOnelabChoices("GmshMerge/InputFiles"));
  std::cout << "Simulation completed successfully..." << std::endl;
  PostArray(Msg::GetOnelabChoices("PostArray"));
}
