#include "OnelabClients.h"

onelab::server *onelab::server::_server = 0;
onelab::remoteNetworkClient *Msg::loader = 0;

void PrintUsage(const char *name){
  printf("\nUsage:       %s [-v int -m int] metamodel\n", name);
  printf("Options are: m      model number (default=0)\n");
  printf("             v      verbosity level (default=0)\n");
  exit(1);
}

int main(int argc, char *argv[]){
  //std::string sockName, clientName;
  PromptUser *globalParam = new PromptUser("glob");

  std::string commandLine="",fileName="";
  int modelNumber=0;
  int verbosity=0;

  int i = 1;
  while(i < argc) {
    if(argv[i][0] == '-') {
      if(!strcmp(argv[i] + 1, "v")) {
	i++;
	//globalParam->setVerbosity(atoi(argv[i]));
	verbosity=atoi(argv[i]);
        i++;        
      }
      else if(!strcmp(argv[i] + 1, "m")) {
	i++;
	modelNumber = (atoi(argv[i]));
        i++;        
      }
      else 
	PrintUsage(argv[0]);
    }
    else {
      if (commandLine.empty())
	commandLine.assign(argv[i]);
      else
	fileName.assign(argv[i]);
      i++;
    }
  }

  std::cout << "command==" << commandLine << std::endl;
  std::cout << "fileName==" << fileName << std::endl;

  globalParam->menu(commandLine, fileName, modelNumber); 
}

