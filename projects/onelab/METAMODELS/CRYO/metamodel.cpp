#include <stdlib.h>
#include <string>
#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *loader = 0;

std::string modelName="";
PromptUser *OL = new PromptUser("onelab");
EncapsulatedGmsh *myMesher = new EncapsulatedGmsh("gmsh");
InterfacedElmer *mySolver = new InterfacedElmer("elmer");
std::ofstream resfile("results.txt");


int main(int argc, char *argv[]){
  bool analyzeOnly=false;
  int modelNumber=0;
  std::string sockName = "";

  getOptions(argc, argv, modelNumber, analyzeOnly, modelName, sockName);
  modelName="cryo";

  if (sockName.size())
    loader = new onelab::remoteNetworkClient("loader", sockName);
  Msg::InitializeOnelab("metamodel",""); // _onelabClient = new onelab::localClient("metamodel");

  if (loader)
    std::cout << "ONELAB: " << Msg::Synchronize_Down(loader) << " parameters downloaded" << std::endl;

  if (analyzeOnly)
    analyze();
  else{
    switch (modelNumber){
    case 1:
      OL->setVerbosity(0);
      OL->setNumber("Tcold",77);
      resfile << "With Tcold =" << OL->getNumber("Tcold") << " K " << std::endl;
      resfile << "The optimum time is " << OL->getNumber("tmin") << " sec." << std::endl;
      compute();
      break;
    case 2: //Computation of the Tcold vs tmin characteristic with a larger time step
      OL->setNumber("NumStep",20); // overrides default value defined in .sif_onelab file
      OL->setNumber("TimeStep",0.1); // idem
      for (int temp=60; temp<120; temp+=10){
	OL->setNumber("Tcold",temp);
	compute();
	resfile << OL->stateToChar() << std::endl; // record -> db
      }
      break;
    case 0:
    default:
      compute();
      break;
    }
  }

  if (loader){
    std::cout << "ONELAB: " << Msg::Synchronize_Up(loader) << " parameters uploaded" << std::endl;
    delete loader;
  }

  Msg::FinalizeOnelab();
}


int analyze(){
  checkIfPresent(modelName+".geo");
  myMesher->analyze("",modelName);
  checkIfPresent(modelName+".sif_onelab");
  mySolver->analyze("",modelName);
  std::cout << OL->showParamSpace() << std::endl;
  return 1;
}

int compute(){

  newStep();
  analyze();

  myMesher->run("-2", modelName);
  checkIfModified(modelName+".msh");
  OL->setString("Gmsh/MshFileName", modelName+".msh");

  std::string cmd="ElmerGrid 14 2 " + modelName + ".msh -out " + OL->getString("elmer/MshDirName");
  systemCall(cmd); // no server access needed
  checkIfModified(OL->getString("elmer/MshDirName")+"/mesh.header");

  mySolver->convert(modelName);
  mySolver->run("",modelName);
  checkIfModified(modelName+".sif");
  checkIfModified("solution.pos");
  checkIfModified("tempevol.txt");

  std::cout << "Simulation completed successfully" << std::endl;

  GmshDisplay(loader,modelName,"solution.pos");
  
  // Post-processing
  cmd="gmsh - solution.pos script.pos.opt"; // compute objective function 
  systemCall(cmd); 

  std::vector<double> v1=extract_column(8,read_array("f1.txt",' '));
  std::vector<double> v2=extract_column(8,read_array("f2.txt",' '));
  std::vector<double> time=extract_column(2,read_array("f1.txt",' '));
  double min=1e200, temp;
  int imin=0;
  
  for (int i=0; i<v1.size(); i++){
    temp = v1[i]+v2[i];
    if (temp<min){ imin=i; min=temp;}
  }
  OL->setNumber("fmin",min,"Value of the objective function");

  if (OL->getNumber("elmer/TRANSIENT")){
    OL->setNumber("tmin",time[imin],"Time at which the objective function is minimum");
    std::cout << std::endl << "Simulation result: tmin=" << time[imin] << std::endl;

    double teval=2.; // evaluation at time=2s
    array data = read_array("tempevol.txt",' '); // temperature over time at probe point
    double T2s = find_in_array((int)(teval/OL->getNumber("elmer/TimeStep")),2,data); // temperature after 2s
  
    OL->setNumber("T2s",T2s,"temperature at (Xloc, Yloc) after 2 sec");
    std::cout << "Simulation result: T2s="  << T2s << std::endl;
  }


  std::cout << OL->showParamSpace() << std::endl;

  return 1;
}
