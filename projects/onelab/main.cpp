#include <stdlib.h>
#include <string>
#include "StringUtils.h"
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *Msg::loader = 0;

int main(int argc, char *argv[]){
  std::string action="", commandLine="", caseName="", clientName="", sockName = "";
  int modelNumber=0;

  getOptions(argc, argv, action, commandLine, caseName, clientName, sockName, modelNumber);
  
  // Msg::_onelabclient is a onelab:LocalClient independent of MetaModel
  Msg::InitializeOnelab("metamodel","");

  if (sockName.size()){
    Msg::loader = new onelab::remoteNetworkClient(clientName, sockName);
    Msg::hasGmsh = clientName.compare("loadedMetaModel");
  }
  else
    Msg::hasGmsh=false;

  Msg::SetOnelabNumber("HasGmsh",Msg::hasGmsh,false);

  if(Msg::loader)
    std::cout << "ONELAB: " << Msg::Synchronize_Down() << " parameters downloaded" << std::endl;

  std::string fileName="", workingDir="";
  if(!caseName.compare("untitled")){
    // no casename was given in calling the metamodel 
    // this means it was called by a loader
    // obtain the caseName from the loader
    caseName= Msg::GetOnelabString("Gmsh/MshFileName");
    if(caseName.empty()){
      caseName= Msg::GetOnelabString("loadedMetaModel/InputFiles");
    }
  }

  std::cout << "caseName=<" << caseName << ">" << std::endl;

  if(caseName.size()){
     fileName.assign(SplitFileName(caseName)[1]);
     workingDir.assign(SplitFileName(caseName)[0]);
     if(workingDir.empty()) workingDir.assign("."+dirSep);
  }
  else
    Msg::Fatal("No valid input filename.");

  Msg::Info("Filename <%s> and working dir <%s>",fileName.c_str(),
	    workingDir.c_str());
  Msg::SetOnelabString("Arguments/FileName",fileName,false);
  Msg::SetOnelabString("Arguments/WorkingDir",workingDir,false);

  MetaModel *myModel = new MetaModel(commandLine, workingDir, clientName,
				     fileName, modelNumber);
  if(!myModel->checkCommandLines()) 
    //if all clients have valid commandlines and are initialized
    action.assign("exit");
  if(Msg::loader && !Msg::GetOnelabNumber(clientName + "/Initialized"))
    action.assign("initialize");

  // std::cout << "checkcmd:" << Msg::GetOnelabString(clientName+"/9CheckCommand") << std::endl;
  // std::cout << "     action:" << action << std::endl;
  // std::cout << "initialized:" << Msg::GetOnelabNumber(clientName + "/Initialized") << std::endl;
  // std::cout << "    hasGmsh:" << Msg::hasGmsh << std::endl;

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

  Msg::Info("Leave metamodel");
}

