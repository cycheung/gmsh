#include <iostream>
#include <vector>
#include <string>
#include "onelab.h"

// 1) Compile this solver: olmake.sh
//
// 2) Add it to Gmsh:
//   - launch Gmsh and open Tools->OneLab 
//   - In the gear menu, select "Add new client"
//   - Enter "My Solver" as client name, then choose the exe in the dialog

int main(int argc, char **argv) 
{
  onelab::remoteNetworkClient *client = 0;

  std::cout << "Running <remote>\n";

  for(int i = 0; i < argc; i++){
    if(std::string(argv[i]) == "-onelab" && i < argc - 2)
      client = new onelab::remoteNetworkClient(argv[i+1],argv[i+2]);
  }

  if(!client){
    printf("usage: %s -onelab clientname socket\n", argv[0]);
    exit(1);
  }
  std::cout << "Remote: I have a client\n";

  std::vector<onelab::string> strings;
  client->get(strings,client->getName()+"/9CheckCommand"); // am I initialized
  if(strings.size()){
    std::cout << "I am initialized: CheckCommand=<" << strings[0].getValue() << ">" << std::endl;
  }
  else{ // initialize
    // send a value to the server
    onelab::string s(client->getName()+"/9CheckCommand","-check");
    client->set(s);
    onelab::number o(client->getName()+"/Initialized",1);
    client->set(o);
    delete client;
    return 0;
  }

  onelab::number o("alpha",123456);
  client->set(o);
  onelab::number oo("beta",654321);
  client->set(oo);

  delete client;
  return 0;
}
