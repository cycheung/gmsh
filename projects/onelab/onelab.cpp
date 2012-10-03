#include <stdlib.h>
#include <string>
#include "StringUtils.h"
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *OLMsg::loader = 0;

int main(int argc, char *argv[]){
  std::string commandLine="", caseName="", clientName="", sockName="";
  parseMode todo;

  getOptions(argc, argv, todo, commandLine, caseName, clientName, sockName);
  
  // OLMsg::_onelabclient is a onelab:LocalClient independent of MetaModel
  OLMsg::InitializeOnelab("meta","");

  if (sockName.size()){
    OLMsg::loader = new onelab::remoteNetworkClient(clientName, sockName);
    OLMsg::hasGmsh = clientName.compare("loadedMetaModel");
  }
  else
    OLMsg::hasGmsh=false;

  OLMsg::SetOnelabNumber("HasGmsh",OLMsg::hasGmsh,false);

  if(OLMsg::loader)
    std::cout << "ONELAB: " << OLMsg::Synchronize_Down() 
	      << " parameters downloaded" << std::endl;

  std::string fileName="", workingDir="";
  if(!caseName.compare("untitled")){
    // no casename was given in calling the metamodel 
    // this means it was called by a loader
    // obtain the caseName from the loader
    caseName= OLMsg::GetOnelabString("Gmsh/MshFileName");
    if(caseName.empty()){
      caseName= OLMsg::GetOnelabString("loadedMetaModel/CaseName");
    }
  }

  if(caseName.size()){
     fileName.assign(SplitFileName(caseName)[1]);
     workingDir.assign(SplitFileName(caseName)[0]);
     //if(workingDir.empty()) workingDir.assign("."+dirSep);
  }
  else
    OLMsg::Fatal("No valid input filename.");

  OLMsg::Info("Filename <%s> and working dir <%s>",fileName.c_str(),
	    workingDir.c_str());
  OLMsg::SetOnelabString("Arguments/FileName",fileName,false);
  OLMsg::SetOnelabString("Arguments/WorkingDir",workingDir,false);

  MetaModel *myModel = new MetaModel(commandLine, workingDir, clientName,
				     fileName);
  myModel->setTodo(todo);

  if(OLMsg::GetOnelabNumber("LOGFILES")){
    freopen("stdout.txt","w",stdout);
    freopen("stderr.txt","w",stderr);
  }

  //if not all clients have valid commandlines -> exit metamodel
  //commandlines will be entered by the user interactively
  if(!myModel->checkCommandLines()) myModel->setTodo(EXIT);
  
  //if the metamodel is not yet initialized by the loader -> initialize
  if(OLMsg::loader && !OLMsg::GetOnelabNumber(clientName + "/Initialized"))
     myModel->setTodo(INITIALIZE);

  if( myModel->isTodo(EXIT)){ 
    // exit metamodel
  }
  else if( myModel->isTodo(INITIALIZE)){
    myModel->initialize(); 
  }
  else if( myModel->isTodo(ANALYZE)){
    myModel->analyze(); 
  }
  else if( myModel->isTodo(COMPUTE)){
    myModel->compute(); 
  }
  else
    OLMsg::Fatal("Main: Unknown Action <%d>", todo);

  if (OLMsg::loader){
    std::cout << "ONELAB: " << OLMsg::Synchronize_Up() 
	      << " parameters uploaded" << std::endl;
    delete OLMsg::loader;
  }
  OLMsg::FinalizeOnelab();
  delete myModel;

  OLMsg::Info("Leave metamodel");
}

