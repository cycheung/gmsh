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
  
  std::vector<std::string> choices;
  std::vector<onelab::string> strings;
  get(strings,"PostArray");
  if(strings.size()){
    choices = strings[0].getChoices();
    int nb=0;
    onelab::number o;
    while( 4*(nb+1) <= choices.size()){
      Msg::Info("PostArray <%s>",choices[4*nb+3].c_str());
      int lin= atof(choices[4*nb+1].c_str());
      int col= atof(choices[4*nb+2].c_str());
      o.setName(choices[4*nb+3]);
      o.setValue(find_in_array(lin,col,read_array(choices[4*nb],' ')));
      set(o);
      nb++;
    }
  }
}
