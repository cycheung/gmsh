#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

//PromptUser *OL = new PromptUser("onelab");
void MetaModel::analyze(){
  Msg::Info("\nMetamodel: now ANALYSING");
  simpleCheck();
  //std::cout << OL->showParamSpace() << std::endl;
  //std::cout << OL->showClientStatus() << std::endl;
}

void MetaModel::compute(){
  Msg::Info("\nMetamodel: now COMPUTING");
  simpleCompute();

  std::vector<std::string> choices;
  if(Msg::GetOnelabChoices("GmshMerge/InputFiles",choices))
    GmshDisplay(Msg::loader,Msg::GetOnelabString("Arguments/FileName"),choices);
  if(Msg::GetOnelabChoices("PostArray",choices))
    PostArray(choices);
}
