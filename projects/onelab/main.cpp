#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *Msg::loader = 0;

// main() commun à tous les métamodèles: elmer file, cryosurgery, etc...
// la classe localClientMetaModel fait l'interface entre les loaders (gmsh, console)
// et les clients solvers

int main(int argc, char *argv[]){
  std::string action="", commandLine="",  fileName="", clientName="", sockName = "";
  int modelNumber=0;

  getOptions(argc, argv, action, commandLine, fileName, clientName, sockName, modelNumber);
  
  // Msg::_onelabclient is a onelab:LocalClient independent of MetaModel

  Msg::InitializeOnelab("metamodel",""); 
  if (sockName.size())
    Msg::loader = new onelab::remoteNetworkClient(clientName, sockName);

  if(Msg::loader)
    std::cout << "ONELAB: " << Msg::Synchronize_Down() << " parameters downloaded" << std::endl;

  if(!fileName.compare("untitled")){ 
    // no filename was given in calling the metamodel => initialize only
    action.assign("initialize");
    // guess fileName from Gmsh
    std::string name= Msg::GetOnelabString("Gmsh/MshFileName");
    if(name.size()){
      fileName.assign(name.substr(0,name.find_last_of(".")));  // remove extension
      Msg::Info("Guessed filename <%s> from OL variable <Gmsh/MshFileName>",fileName.c_str());
    }
    else
      Msg::Info("No valid input filename found");
  }
  Msg::SetOnelabString("Arguments/FileName",fileName);

  MetaModel *myModel = new MetaModel(commandLine, clientName, fileName, modelNumber);

  newStep();
  if(!action.compare("initialize")){
    if(Msg::loader) myModel->initialize();
  }
  else if(!action.compare("check")){
    myModel->analyze();
  }
  else if(!action.compare("compute")){
    myModel->analyze();
    myModel->compute();
  }
  else
    Msg::Fatal("Main: Unknown Action <%s>", action.c_str());

  if (Msg::loader){
    std::cout << "ONELAB: " << Msg::Synchronize_Up() << " parameters uploaded" << std::endl;
    delete Msg::loader;
  }
  Msg::FinalizeOnelab();
  delete myModel;

  std::cout << "ONELAB: leave metamodel." << std::endl;
}

