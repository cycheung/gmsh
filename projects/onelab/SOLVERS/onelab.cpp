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
  Msg::Info("Simulation completed successfully...");
  GmshDisplay(Msg::loader,Msg::GetOnelabString("Arguments/FileName"),Msg::GetOnelabChoices("GmshMerge/InputFiles"));
  PostArray(Msg::GetOnelabChoices("PostArray"));
  Msg::Info("Post-processing completed successfully...");
}
