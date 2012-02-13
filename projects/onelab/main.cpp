#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *Msg::loader = 0;
//std::set<std::string, ShortNameLessThan> localSolverClient::_parameters;

// main() commun à tous les métamodèles

int main(int argc, char *argv[]){
  std::string action="", commandLine="",  fileName="", clientName="", sockName = "";
  int modelNumber=0;

  getOptions(argc, argv, action, commandLine, fileName, clientName, sockName, modelNumber);
  
  // Msg::_onelabclient is a onelab:LocalClient independent of MetaModel
  Msg::InitializeOnelab("metamodel","");

  if (sockName.size()){
    Msg::loader = new onelab::remoteNetworkClient(clientName, sockName);
    //Msg::hasGmsh = loader->getName().compare("loadedMetaModel");
    Msg::hasGmsh = clientName.compare("loadedMetaModel");
  }

  if(Msg::loader)
    std::cout << "ONELAB: " << Msg::Synchronize_Down() << " parameters downloaded" << std::endl;

  if(!fileName.compare("untitled")){ 
    // no filename was given in calling the metamodel 
    // try to guess 'Arguments/FileName' from arguments of the loader
    std::string name= Msg::GetOnelabString("Gmsh/MshFileName");
    if(name.size()){
      fileName.assign(name.substr(0,name.find_last_of(".")));  // remove extension
      Msg::Info("Guessed filename <%s> from OL variable <Gmsh/MshFileName>",fileName.c_str());
    }
    else{
      name= Msg::GetOnelabString("loadedMetaModel/InputFiles");
      if(name.size()){
	fileName.assign(name.substr(0,name.find_last_of(".")));  // remove extension
	Msg::Info("Guessed filename <%s> from OL variable <loadedMetaModel/InputFiles>",fileName.c_str());
      }
      else{
	Msg::Info("No valid input filename found, action set to <initialize>");
      }
    }
  }
  Msg::SetOnelabString("Arguments/FileName",fileName,false);

  //std::cout << "fileName:" << fileName << std::endl;

  MetaModel *myModel = new MetaModel(commandLine, clientName, fileName, modelNumber);
  if(!myModel->checkPathes()) //true if all clients have path and are initialized
    action.assign("exit");
  if(Msg::loader && !Msg::GetOnelabNumber(clientName + "/Initialized"))
    action.assign("initialize");

  std::cout << "   checkcmd:" << Msg::GetOnelabString(clientName+"/9CheckCommand") << std::endl;
  std::cout << "     action:" << action << std::endl;
  std::cout << "initialized:" << Msg::GetOnelabNumber(clientName + "/Initialized") << std::endl;
  std::cout << "    hasGmsh:" << Msg::hasGmsh << std::endl;

  //newStep();
  if(!action.compare("exit")){ // exit metamodel
  } 
  else if(!action.compare("initialize")){
    if(Msg::loader) myModel->initialize(); // initializes MetaModel 
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

