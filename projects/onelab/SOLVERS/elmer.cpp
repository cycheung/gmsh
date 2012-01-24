#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

PromptUser *OL = new PromptUser("onelab");
EncapsulatedClient *myMesher = new EncapsulatedClient("Gmsh","gmsh",".geo");
InterfacedClient *mySolver = new InterfacedClient("Elmer","ElmerSolver",".sif");

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

  std::string cmd="ElmerGrid 14 2 " + genericNameFromArgs + ".msh -out " + OL->getString(mySolver->getName() + "/MshDirName");
  systemCall(cmd); // no server access needed
  checkIfModified(OL->getString(mySolver->getName() + "/MshDirName")+"/mesh.header");

  mySolver->setLineOptions("");
  mySolver->convert();
  mySolver->compute();

  GmshDisplay(Msg::loader,genericNameFromArgs,Msg::GetOnelabChoices(mySolver->getName() + "/9OutputFiles"));

  std::cout << "Simulation completed successfully" << std::endl;

  array data = read_array("tempevol.txt",' '); // temperature over time at probe point
  double Tpost = find_in_array(1,2,data); 
  //OL->setNumber("Tpost",Tpost,"temperature at (Xloc, Yloc)");
  OL->AddNumberChoice("Tpost",Tpost);
  std::cout << "Simulation result: Tpost="  << Tpost << std::endl;
}
