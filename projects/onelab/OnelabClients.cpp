#include "OnelabMessage.h"
#include "OnelabClients.h"
#include "StringUtils.h"
#include <algorithm>


class onelabMetaModelServer : public GmshServer{
 private:
  localNetworkSolverClient *_client;
 public:
  onelabMetaModelServer(localNetworkSolverClient *client)
    : GmshServer(), _client(client) {}
  ~onelabMetaModelServer(){}
  int NonBlockingSystemCall(const char *command)
  {
#if defined(WIN32)
    STARTUPINFO suInfo;
    PROCESS_INFORMATION prInfo;
    memset(&suInfo, 0, sizeof(suInfo));
    suInfo.cb = sizeof(suInfo);
    std::string cmd(command);
    Msg::Info("Calling <%s>", cmd.c_str());
    // DETACHED_PROCESS removes the console (useful if the program to launch is
    // a console-mode exe)
    CreateProcess(NULL,(char *)cmd.c_str(), NULL, NULL, FALSE,
		  NORMAL_PRIORITY_CLASS|DETACHED_PROCESS, NULL, NULL,
		  &suInfo, &prInfo);
    return 0;
#else
    if(!system(NULL)) {
      Msg::Error("Could not find /bin/sh: aborting system call");
      return 1;
    }
    std::string cmd(command);
    int pos;
    if((pos=cmd.find("incomp_ssh ")) != std::string::npos){
      cmd.assign(cmd.substr(pos+7));  // remove "incomp_"  
      cmd.append(" & '");
    }
    else 
      cmd.append(" & ");

    Msg::Info("Calling <%s>", cmd.c_str());
    return system(cmd.c_str());
#endif
  }// non blocking 
  int NonBlockingWait(int socket, double waitint, double timeout)
  {
    double start = GetTimeInSeconds();
    while(1){
      if(timeout > 0 && GetTimeInSeconds() - start > timeout)
        return 2; // timout
      // if(_client->getPid() < 0 || (_client->getCommandLine().empty() &&
      //                              !CTX::instance()->solver.listen))
      if(_client->getPid() < 0 || (_client->getCommandLine().empty()))
        return 1; // process has been killed or we stopped listening
      // check if there is data (call select with a zero timeout to
      // return immediately, i.e., do polling)
      int ret = Select(0, 0, socket);
      if(ret == 0){ // nothing available
        // if asked, refresh the onelab GUI
        std::vector<onelab::string> ps;
        onelab::server::instance()->get(ps, "Gmsh/Action");
        if(ps.size() && ps[0].getValue() == "refresh"){
          ps[0].setVisible(false);
          ps[0].setValue("");
          onelab::server::instance()->set(ps[0]);
          //onelab_cb(0, (void*)"refresh");
        }
        // wait at most waitint seconds and respond to FLTK events
        //FlGui::instance()->wait(waitint);
      }
      else if(ret > 0){
        return 0; // data is there!
      }
      else{
        // an error happened
        _client->setPid(-1);
        _client->setGmshServer(0);
        return 1;
      }
    }
  }
};

std::string localNetworkSolverClient::buildCommandLine(){
  std::string command;
  command.assign("");
  if(!getWorkingDir().empty()) 
    command.append("cd " + getWorkingDir() + cmdSep);
  command.append(FixWindowsPath(getCommandLine()));
  if(command.size()){
    std::vector<onelab::string> ps;
    get(ps, getName() + "/Action");
    std::string action = (ps.empty() ? "" : ps[0].getValue());
    get(ps, getName() + "/9CheckCommand");
    std::string checkCommand = (ps.empty() ? "" : ps[0].getValue());
    get(ps, getName() + "/9ComputeCommand");
    std::string computeCommand = (ps.empty() ? "" : ps[0].getValue());

    if(action == "initialize"){
      command.append(" " + getSocketSwitch() + " " + getName() + " %s"); 
    }
    else if(action == "check") {
      command.append(" " + getString("Arguments") + " " + checkCommand) ;
      command.append(" " + getSocketSwitch() + " " + getName() + " %s"); 
    }
    else if(action == "compute"){
      command.append(" " + getString("Arguments") + " " + computeCommand);
      command.append(" " + getSocketSwitch() + " " + getName() + " %s"); 
      command.append(" " + getString("Redirect"));
      //does not quite work. Maybe using fork instead of system would help...
    }
    else
      Msg::Fatal("localNetworkSolverClient::buildCommandLine: Unknown Action <%s> %s", action.c_str(), getName().c_str());
  }
  return command;
}

bool localNetworkSolverClient::run()
{
 new_connection:
  _pid = 0;
  _gmshServer = 0;

  onelabMetaModelServer *server = new onelabMetaModelServer(this);
 
#if defined(WIN32)
  std::string socketName = ":";
#else
  std::string socketName;
  if(getRemote())
    socketName = ":";
  else
    socketName = getUserHomedir() + ".gmshsock";
#endif

  std::string sockname;
  std::ostringstream tmp;
  if(!strstr(socketName.c_str(), ":")){
    // Unix socket
    tmp << socketName << getId();
    sockname = FixWindowsPath(tmp.str());
  }
  else{
    // TCP/IP socket
    if(socketName.size() && socketName[0] == ':')
      // prepend hostname if only the port number is given
      tmp << GetHostName(); 
    tmp << socketName ;
    sockname = tmp.str();
  }

  std::string command = buildCommandLine();
  // std::cout << "sockname=<" << sockname << ">" << std::endl;
  // std::cout << "command=<" << command << ">" << std::endl;

  int sock;
  try{
    sock = server->Start(command.c_str(), sockname.c_str(), 10);
  }
  catch(const char *err){
    Msg::Error("%s (on socket '%s')", err, sockname.c_str());
    sock = -1;
  }

  if(sock < 0){
    server->Shutdown();
    delete server;
    return false;
  }

  Msg::StatusBar(2, true, "Now running client <%s>", _name.c_str());

  while(1) {
    if(_pid < 0) break;
    
    int stop = server->NonBlockingWait(sock, 0.1, 0.);
    if(stop || _pid < 0) {
      Msg::Info("Stop=%d _pid=%d",stop, _pid);
      break;
    }
    int type, length, swap;
    if(!server->ReceiveHeader(&type, &length, &swap)){
      Msg::Error("Did not receive message header: stopping server");
      break;
    }
     // else
     //   std::cout << "FHF: Received header=" << type << std::endl;

    std::string message(length, ' ');
    if(!server->ReceiveMessage(length, &message[0])){
      Msg::Error("Did not receive message body: stopping server");
      break;
    }
    // else
    //   std::cout << "FHF: Received message=" << message << std::endl;

    switch (type) {
    case GmshSocket::GMSH_START:
      _pid = atoi(message.c_str());
      _gmshServer = server;
      break;
    case GmshSocket::GMSH_STOP:
      _pid = -1;
      _gmshServer = 0;
      break;
    case GmshSocket::GMSH_PARAMETER:
      {
        std::string version, type, name;
        onelab::parameter::getInfoFromChar(message, version, type, name);
        if(onelab::parameter::version() != version){
          Msg::Error("OneLab version mismatch for %s (server: %s / client: %s)",
	   message.c_str(), onelab::parameter::version().c_str(), version.c_str());
        }
        if(type == "number"){
          onelab::number p; p.fromChar(message); set(p);
        }
        else if(type == "string"){
          onelab::string p; p.fromChar(message); set(p);
        }
	else if(type == "region"){
          onelab::region p; p.fromChar(message); set(p);
        }
        else if(type == "function"){
          onelab::function p; p.fromChar(message); set(p);
        }
        else
          Msg::Fatal("FIXME query not done for this parameter type: <%s>",
		     message.c_str());
      }
      break;
    case GmshSocket::GMSH_PARAMETER_QUERY:
      {
        std::string version, type, name, reply;
        onelab::parameter::getInfoFromChar(message, version, type, name);
        if(onelab::parameter::version() != version){
          Msg::Error("OneLab version mismatch for %s (server: %s / client: %s)",
	   message.c_str(), onelab::parameter::version().c_str(), version.c_str());
        }
        else if(type == "number"){
          std::vector<onelab::number> par; get(par, name);
          if(par.size() == 1) reply = par[0].toChar();
	}
        else if(type == "string"){
          std::vector<onelab::string> par; get(par, name);
          if(par.size() == 1) reply = par[0].toChar();
        }
        else if(type == "region"){
          std::vector<onelab::region> par; get(par, name);
          if(par.size() == 1) par[0].toChar();
        }
        else if(type == "function"){
          std::vector<onelab::function> par; get(par, name);
          if(par.size() == 1) reply = par[0].toChar();
        }
        else
          Msg::Error("Unknown OneLab parameter type in query: %s",
		     type.c_str());

        if(reply.size()){
          server->SendMessage(GmshSocket::GMSH_PARAMETER, 
			      reply.size(), &reply[0]);
        }
        else{
          reply = "OneLab parameter '" + name + "' not found";
          server->SendMessage(GmshSocket::GMSH_INFO, reply.size(), &reply[0]);
        }
      }
      break;
    case GmshSocket::GMSH_PARAM_QUERY_ALL:
      {
        std::string version, type, name, reply;
        std::vector<std::string> replies;
        onelab::parameter::getInfoFromChar(message, version, type, name);
	if(onelab::parameter::version() != version){
          Msg::Error("OneLab version mismatch for %s (server: %s / client: %s)",
           message.c_str(), onelab::parameter::version().c_str(), version.c_str());
        }
	else if(type == "number"){
	  std::vector<onelab::number> numbers; get(numbers);
	  for(std::vector<onelab::number>::iterator it = numbers.begin(); 
	      it != numbers.end(); it++) replies.push_back((*it).toChar());
	}
	else if(type == "string"){
	  std::vector<onelab::string> strings; get(strings);
	  for(std::vector<onelab::string>::iterator it = strings.begin(); 
	      it != strings.end(); it++) replies.push_back((*it).toChar());
	}
        else if(type == "region"){
	  std::vector<onelab::region> regions; get(regions);
	  for(std::vector<onelab::region>::iterator it = regions.begin();
              it != regions.end(); it++) replies.push_back((*it).toChar());
        }
        else if(type == "function"){
	  std::vector<onelab::function> functions; get(functions);
	  for(std::vector<onelab::function>::iterator it = functions.begin();
              it != functions.end(); it++) replies.push_back((*it).toChar());
        }
        else
          Msg::Error("Unknown OneLab parameter type in query: %s",
		     type.c_str());

        for(unsigned int i = 0; i < replies.size(); i++)
          server->SendMessage(GmshSocket::GMSH_PARAM_QUERY_ALL, 
	     replies[i].size(), &replies[i][0]);
        reply = "Sent all OneLab " + type + "s";
        server->SendMessage(GmshSocket::GMSH_PARAM_QUERY_END, 
			    reply.size(), &reply[0]);
      }
      break;
    case GmshSocket::GMSH_PROGRESS:
      Msg::StatusBar(2, false, "%s %s", _name.c_str(), message.c_str());
      break;
    case GmshSocket::GMSH_INFO:
      Msg::Direct("%-8.8s: %s", _name.c_str(), message.c_str());
      break;
    case GmshSocket::GMSH_WARNING:
      Msg::Direct(2, "%-8.8s: %s", _name.c_str(), message.c_str());
      break;
    case GmshSocket::GMSH_ERROR:
      //Msg::Direct(1, "%-8.8s: %s", _name.c_str(), message.c_str());
      Msg::Fatal("%-8.8s: %s", _name.c_str(), message.c_str());
      break;
    case GmshSocket::GMSH_MERGE_FILE:
      Msg::Info("Merge Post-Processing File %s",message.c_str());
      SystemCall("gmsh "+ message);
      //implémentation pour le cas du loader en mode console.
      break;
    default:
      Msg::Warning("Received unknown message type (%d)", type);
      break;
    }
  }

  server->Shutdown();
  delete server;
  Msg::StatusBar(2, true, "Done running '%s'", _name.c_str());
  return true;
}

bool localNetworkSolverClient::kill()
{
  if(_pid > 0) {
    if(KillProcess(_pid)){
      Msg::Info("Killed '%s' (pid %d)", _name.c_str(), _pid);
      _pid = -1;
      return true; 
    }
  }
  _pid = -1;
  return false;
}

// client PROMPTUSER

void PromptUser::setNumber(const std::string paramName, const double val, const std::string &help){
  onelab::number number;
  number.setName(paramName);
  number.setValue(val);
  number.setHelp(help);
  set(number);
}

void PromptUser::addNumberChoice(std::string name, double val){
  std::vector<double> choices;
  std::vector<onelab::number> numbers;
  get(numbers, name);
  if(numbers.size())
    choices = numbers[0].getChoices();
  else{
    numbers.resize(1);
    numbers[0].setName(name);
  }
  numbers[0].setValue(val);
  choices.push_back(val);
  numbers[0].setChoices(choices);
  set(numbers[0]);
}

void PromptUser::addStringChoice(std::string name, std::string str){
  std::vector<std::string> choices;
  std::vector<onelab::string> strings;
  get(strings, name);
  if(strings.size())
    choices = strings[0].getChoices();
  else{
    strings.resize(1);
    strings[0].setName(name);
  }
  strings[0].setValue(str);
  choices.push_back(str);
  strings[0].setChoices(choices);
  set(strings[0]);
}

double PromptUser::getNumber(const std::string paramName){
  std::vector<onelab::number> numbers;
  get(numbers,paramName);
  if (numbers.size())
    return numbers[0].getValue();
  else
    Msg::Fatal("Unknown parameter %s",paramName.c_str());
}
bool PromptUser::existNumber(const std::string paramName){
  std::vector<onelab::number> numbers;
  get(numbers,paramName);
  return numbers.size();
}
void PromptUser::setString(const std::string paramName, const std::string &val, const std::string &help){
  onelab::string string;
  string.setName(paramName);
  string.setValue(val);
  string.setHelp(help);
  set(string);
}
std::string PromptUser::getString(const std::string paramName){
  std::vector<onelab::string> strings;
  get(strings,paramName);
  if (strings.size())
    return strings[0].getValue();
  else
    Msg::Fatal("Unknown parameter %s",paramName.c_str());
}
bool PromptUser::existString(const std::string paramName){
  std::vector<onelab::string> strings;
  get(strings,paramName);
  return strings.size();
}

std::string PromptUser::stateToChar(){
  std::vector<onelab::number> numbers;
  std::ostringstream sstream;
  get(numbers);
  for(std::vector<onelab::number>::iterator it = numbers.begin();
      it != numbers.end(); it++)
    sstream << (*it).getValue() << '\t';
  return sstream.str();
}

std::string PromptUser::showParamSpace(){
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  std::vector<onelab::region> regions;
  std::ostringstream sstream;

  get(numbers);
  sstream << std::endl << itoa(numbers.size()) << " numbers" << std::endl;
  for(std::vector<onelab::number>::iterator it = numbers.begin();
      it != numbers.end(); it++)
    sstream << (*it).getName() << " = " << (*it).getValue() << std::endl;

  get(strings);
  sstream << std::endl << itoa(strings.size()) << " strings" << std::endl;
  for(std::vector<onelab::string>::iterator it = strings.begin();
      it != strings.end(); it++)
    sstream << (*it).getName() << " = " << (*it).getValue() << std::endl;

  get(regions);
  sstream << std::endl << itoa(regions.size()) << " regions" << std::endl;
  for(std::vector<onelab::region>::iterator it = regions.begin();
      it != regions.end(); it++){
    std::set<std::string> region=  (*it).getValue();
    sstream << "<" ;
    for(std::set<std::string>::const_iterator iter = region.begin();
	iter != region.end(); iter++)
      sstream << sanitize(*iter);
    sstream << ">" ;
  }

  return sstream.str();

  // std::string db = "ONELAB parameter space: size=" 
  //   + itoa(onelab::server::instance()->getNumParameters()) + "\n";
  // db.append(onelab::server::instance()->toChar());
  // for(unsigned int i = 0; i < db.size(); i++)
  //   if(db[i] == onelab::parameter::charSep()) db[i] = '|';
  // return db.c_str();

}

std::string PromptUser::showClientStatus(){
  std::ostringstream sstream;
  std::cout << "\nONELAB: Present state of the onelab clients" << std::endl;
  for(onelab::server::citer it = onelab::server::instance()->firstClient();
      it != onelab::server::instance()->lastClient(); it++){
    std::string name= it->second->getName();
    sstream << "<" << onelab::server::instance()->getChanged(name) << "> " << name << std::endl;
  }
  return sstream.str();
}

bool PromptUser::menu(std::string commandLine, std::string workingDir, std::string fileName, int modelNumber) { 
  int choice, counter1=0, counter2=0;
  std::string answ;
  std::string name;
  onelab::number x;
  onelab::string y;
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;

  std::string clientName="loadedMetaModel";
  EncapsulatedClient *loadedSolver = new  EncapsulatedClient(clientName,commandLine,workingDir);
  //setString("Arguments/FileName",fileName);
  addStringChoice(clientName + "/CaseName", workingDir+fileName);

  do {
    std::cout << "\nONELAB: menu" << std::endl ;
    std::cout << " 1- View parameter space\n 2- Set a value\n 3- Analyze\n 4- Compute\n 5- List modified files\n 6- Client status\n 7- Quit metamodel" << std::endl;
    choice=0;
    std::string mystr;
    while( (choice<1 || choice>7) && ++counter1<10 ) {
      std::cout << "\nONELAB: your choice? "; 
      std::getline (std::cin, mystr);
      std::stringstream myStream(mystr);
      if (myStream >> choice) break;
      std::cout << "Invalid choice" << std::endl;
    }
    std::cout << "Your choice is <" << choice << ">" << std::endl;

    if (choice==1){
      std::cout << showParamSpace();
      choice=0;
    }
    else if (choice==2){
      std::cout << "ONELAB: Variable name? "; std::cin >> name;
      get(numbers,name);
      if (numbers.size()) {
	float fval;
	std::cout << "ONELAB: Value? "; std::cin >> fval;
	numbers[0].setValue(fval);
	bool allowed = set(numbers[0]);
      }
      else{
	get(strings,name);
	if (strings.size()) {
	  std::string sval;
	  std::cout << "ONELAB: Value? "; std::cin >> sval;
	  strings[0].setValue(sval);
	  bool allowed = set(strings[0]);
	  std::cout << "ONELAB: Ok you are allowed to modify this!" << std::endl;
	}
	else
	  std::cout << "ONELAB: The variable " << name << " is not defined" << std::endl;
      }
      choice=0;
    }
    else if (choice==3){
      loadedSolver->analyze();
      choice=0;
    }
    else if (choice==4){
      loadedSolver->compute();
      choice=0;
    }
    else if (choice==5){
      std::ifstream infile("onelab.modified");
      std::string buff;
      if (infile.is_open()){
	while ( infile.good() ) {
	  getline (infile,buff);
	  std::cout << buff << std::endl;
	}
      }
      choice=0;
    }
    else if (choice==6){
      std::cout << "\nONELAB: Present state of the onelab clients\n" << std::endl;
      for(onelab::server::citer it = onelab::server::instance()->firstClient();
	  it != onelab::server::instance()->lastClient(); it++){
	std::string name= it->second->getName();
	std::cout << name << "<" << onelab::server::instance()->getChanged(name) << ">" << std::endl;
      }
      choice=0;
    }
    else if (choice==7)
      exit(1);
    else
      choice=0;
  } while(!choice && ++counter2<20);
}

// client LOCALSOLVERCLIENT

std::string localSolverClient::toChar(){
  std::ostringstream sstream;
  if(getCommandLine().size()){
    sstream << getName() << "." << "commandLine(" 
	    << getCommandLine() << ");\n";
  }
  return sstream.str();
}

bool localSolverClient::checkCommandLine(){
  Msg::Info("Check command line for <%s>",getName().c_str());
  if(getCommandLine().empty()){
    // look if one has a commandLine on server
    std::string commandLine = getString("CommandLine");
    if(!commandLine.empty()){
      setCommandLine(commandLine);
    }
  }

  if(!isActive()) return true;

  if(!getCommandLine().empty()){
    Msg::SetOnelabString(getName() + "/Action","initialize",false);
    run();
  }
  else{
    if(Msg::hasGmsh) {
      // exits metamodel and restores control to the onelab window
      Msg::Error("The command line of client <%s> is undefined.", getName().c_str());
      std::cout << "\n\nBrowse for executable in the ONELAB window.\n\n" << std::endl;
      return false;
    }
    else{ // asks the user in console mode
      std::cout << "\nONELAB:Enter the command line (with path) of the executable file of <" << getName() << ">" << std::endl;
      std::string cmdl;
      std::getline (std::cin,cmdl);
      setCommandLine(cmdl);
      return cmdl.size();
    }
  }
  return true;
}

const std::string localSolverClient::getString(const std::string what){
  std::vector<onelab::string> strings;
  std::string paramName(getName()+"/"+what);
  get(strings,paramName);
  if(strings.size())
    return strings[0].getValue();
  else
    return "";
}

const bool localSolverClient::getList(const std::string type, std::vector<std::string> &choices){
  std::vector<onelab::string> strings;
  get(strings,getName()+"/"+type);
  if(strings.size()){
    choices= strings[0].getChoices();
    return true;
  }
  else
    return false;
}

bool localSolverClient::buildRmCommand(std::string &cmd){
  std::vector<std::string> choices;

  cmd.assign("");
  if(!getWorkingDir().empty()) 
    cmd.append("cd " + getWorkingDir() + cmdSep);

  if(getList("OutputFiles",choices)){
#if defined(WIN32)
    cmd.append("del ");
#else
    cmd.append("rm -rf ");
#endif
    for(unsigned int i = 0; i < choices.size(); i++)
      cmd.append(choices[i]+" ");
    return true;
  }
  else
    return false;
}

void localSolverClient::PostArray(std::vector<std::string> choices)
{
  int nb=0;
  onelab::number o;
  while( 4*(nb+1) <= choices.size()){
    //std::cout << "Nb Choices" << choices.size() << std::endl;
    int lin= atof(choices[4*nb+1].c_str());
    int col= atof(choices[4*nb+2].c_str());
    std::string fileName = getWorkingDir()+choices[4*nb];
    double val=find_in_array(lin,col,read_array(fileName,' '));
    Msg::AddOnelabNumberChoice(choices[4*nb+3],val);
    Msg::Info("Upload parameter <%s>=%e from file <%s>", choices[4*nb+3].c_str(),val,fileName.c_str());
    nb++;
  }
}

void localSolverClient::GmshMerge(std::vector<std::string> choices)
{
  if(!Msg::hasGmsh || choices.empty()) return;
  for(unsigned int i = 0; i < choices.size(); i++){
    if(!choices[i].empty()){
      std::string fileName=getWorkingDir()+choices[i];
      Msg::loader->sendMergeFileRequest(fileName);
      Msg::Info("Send merge request <%s>",fileName.c_str());
    }
  }
}

// client METAMODEL

//Metamodel::analyze and Metamodel::compute are defined in the file SOLVERS/onelab.cpp

void MetaModel::saveCommandLines(const std::string fileName){ 
  //save client command lines
  std::string fileNameSave = getWorkingDir()+fileName+onelabExtension+".save";
  std::ofstream outfile(fileNameSave.c_str()); 

  if (outfile.is_open())
    for(citer it = _clients.begin(); it != _clients.end(); it++)
      outfile << (*it)->toChar();
  else
    Msg::Fatal("The file <%s> cannot be opened",fileNameSave.c_str());
  outfile.close();
}

bool MetaModel::checkCommandLines(){
  bool allDefined=true;
  for(citer it = _clients.begin(); it != _clients.end(); it++){
    allDefined = allDefined && (*it)->checkCommandLine();
  }
  saveCommandLines(genericNameFromArgs);
  return allDefined;
}

void MetaModel::initialize()
{
  Msg::Info("Initialize Metamodel by the loader");
  Msg::SetOnelabString(clientName + "/9CheckCommand","-a",false);
  Msg::SetOnelabNumber(clientName + "/UseCommandLine",1,false);
  Msg::SetOnelabNumber(clientName + "/Initialized",1,false);
}

void MetaModel::registerClient(const std::string &name, const std::string &type, const std::string &cmdl, const std::string &host, const std::string &rdir) {
  localSolverClient *c;
  // Clients are assigned by default the same working dir as the MetaModel
  // i.e. the working dir from args
  if(host.empty() || rdir.empty()){ //local client
    if(!type.compare(0,6,"interf"))
      c= new InterfacedClient(name,cmdl,getWorkingDir());
    else if(!type.compare(0,6,"encaps"))
      c= new EncapsulatedClient(name,cmdl,getWorkingDir());
    else 
      Msg::Fatal("Unknown client type", type.c_str());
  }
  else{ // remote client
    if(!type.compare(0,6,"interf"))
      c= new RemoteInterfacedClient(name,cmdl,getWorkingDir(),host,rdir);
    else if(!type.compare(0,6,"encaps"))
      c= new RemoteEncapsulatedClient(name,cmdl,getWorkingDir(),host,rdir);
    else 
      Msg::Fatal("Unknown remote client type", type.c_str());
  }
  _clients.push_back(c);
}

void MetaModel::analyze()
{
  Msg::Info("\nMetamodel: now ANALYSING");
  for(citer it = _clients.begin(); it != _clients.end(); it++){
    //if((*it)->isActive()){ // also inactive clients are checked
      Msg::SetOnelabString((*it)->getName() + "/Action","check",false);
      Msg::Info("%d <%s>", 
		onelab::server::instance()->getChanged((*it)->getName()),
		(*it)->getName().c_str());
      (*it)->analyze();

      //some clients must be run at the check phase
      if((*it)->getActive()==2) (*it)->compute();
  }
}

void MetaModel::compute()
{
  Msg::Info("\nMetamodel: now COMPUTING");
  for(citer it = _clients.begin(); it != _clients.end(); it++){
    if(Msg::GetOnelabString("MetaModel/STATUS").compare("STOP")){
      if((*it)->isActive()){
	Msg::SetOnelabString((*it)->getName() + "/Action","compute",false);
	// if(Msg::GetOnelabNumber("LOGFILES")){
	//   freopen((*it)->getName().append(".log").c_str(),"w",stdout);
	//   freopen((*it)->getName().append(".err").c_str(),"w",stderr);
	// }
	(*it)->compute();
	Msg::SetOnelabString((*it)->getName() + "/Action","alldone",false);
      }
    }
  }
  // resets all parameters to 'unchanged'
  std::cout << "FHF sets all param unchanged" << std::endl;
  onelab::server::instance()->setChanged(false);
}

void MetaModel::PostArray(std::vector<std::string> choices)
{
  int nb=0;
  onelab::number o;
  while( 4*(nb+1) <= choices.size()){
    int lin= atof(choices[4*nb+1].c_str());
    int col= atof(choices[4*nb+2].c_str());
    std::string filename = Msg::GetOnelabString("Arguments/WorkingDir")+choices[4*nb];
    double val=find_in_array(lin,col,read_array(filename,' '));
    Msg::AddOnelabNumberChoice(choices[4*nb+3],val);
    Msg::Info("PostArray <%s>=%e",choices[4*nb+3].c_str(),val);
    nb++;
  }
}

// INTERFACED client

void InterfacedClient::analyze() {
  int pos;
  std::vector<std::string> choices;
  Msg::SetOnelabString(getName() + "/Action","check",false);
  getList("InputFiles", choices);
  for(unsigned int i = 0; i < choices.size(); i++){
    if((pos=choices[i].find(onelabExtension)) != std::string::npos){
      // if .ol file
      checkIfPresentLocal(choices[i]);
      parse_onefile(choices[i]);
    }
  }
  convert(); // only for debug
}

void InterfacedClient::convert() {
  int pos;
  std::vector<std::string> choices;
  getList("InputFiles", choices);
  for(unsigned int i = 0; i < choices.size(); i++){
    if((pos=choices[i].find(onelabExtension)) != std::string::npos){
      checkIfPresentLocal(choices[i]);
      // remove .ol extension
      std::string ofilename = getWorkingDir() + choices[i].substr(0,pos);
      std::ofstream outfile(ofilename.c_str());
      if (outfile.is_open())
	convert_onefile(choices[i],outfile);
      else
	Msg::Fatal("The file <%s> cannot be opened",ofilename.c_str());
      outfile.close();
    }
  }
}

void InterfacedClient::compute(){
  std::string cmd;
  std::vector<std::string> choices;
  analyze();
  convert(); 
  Msg::SetOnelabString(getName() + "/Action","compute",false);
  if(getList("InputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++){
     //remove .ol ext
     checkIfPresentLocal(choices[i].substr(0,choices[i].find(onelabExtension)));
    }
  }

  if(getActive())
    if(buildRmCommand(cmd))
      SystemCall(cmd,true); //blocking

  cmd.assign("");
  if(!getWorkingDir().empty())
    cmd.append("cd " + getWorkingDir() + cmdSep);
  cmd.append(FixWindowsPath(getCommandLine() + " "));
  cmd.append(getString("Arguments"));

  SystemCall(cmd.c_str(),true); //blocking

  if(getList("OutputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++){
      checkIfPresentLocal(choices[i]);
    }
  }

  if(getList("PostArray",choices))
    PostArray(choices);

  if(getList("Merge",choices))
    GmshMerge(choices);

  Msg::Info("Client %s completed",_name.c_str());
}

// ENCAPSULATED Client

void EncapsulatedClient::analyze() {
  Msg::SetOnelabString(getName() + "/Action","check",false);
  run();
}

void EncapsulatedClient::compute() {
  std::string cmd;
  std::vector<std::string> choices;

  Msg::SetOnelabString(getName() + "/Action","compute",false);
  if(getActive())
    if(buildRmCommand(cmd))
      SystemCall(cmd,true);

  if(getList("InputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++){
      checkIfPresentLocal(choices[i]);
    }
  }
  run();

  if(getList("OutputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++){
      checkIfPresentLocal(choices[i]);
    }
  }

  if(getList("PostArray",choices))
    PostArray(choices);

  if(getList("Merge",choices))
    GmshMerge(choices);
}

// REMOTE CLIENT

int mySystem(std::string commandLine){
  //std::cout << "mySystem<" << commandLine << ">" << std::endl;
  return SystemCall(commandLine.c_str(), true);
}

bool remoteClient::checkIfPresentRemote(const std::string &fileName){
  struct stat buf;
  std::string cmd;
  char cbuf [1024];
  FILE *fp;

  cmd.assign("ssh "+_remoteHost+" 'cd "+_remoteDir+"; ls "+fileName+" 2>/dev/null'");
  //std::cout << "check remote<" << cmd << ">" << std::endl;
  fp = popen(cmd.c_str(), "r");
  if(fgets(cbuf, 1024, fp) == NULL){
    Msg::Fatal("The file <%s> is not present", fileName.c_str());
    pclose(fp);
    return false;
  }
  Msg::Info("The file <%s> is present", fileName.c_str());
  pclose(fp);
  return true;
}

bool remoteClient::syncInputFile(const std::string &wdir, const std::string &fileName){
  int pos;
  std::string cmd;
  if((pos=fileName.find(onelabExtension)) != std::string::npos){ 
    // .ol file => local
    std::string trueName = 
      fileName.substr(fileName.find_first_not_of(" "),pos); // remove ext .ol
    std::string fullName = wdir+trueName;
    if(checkIfPresent(fullName)){
      cmd.assign("rsync -e ssh -auv "+fullName+" "+_remoteHost+":"+_remoteDir+"/"+trueName);
      return mySystem(cmd);
    }
    else{
      Msg::Fatal("The input file <%s> is not present", fullName.c_str());
      return false;
    }
  }
  else { // not a .ol file, apply the . rule
    if(!fileName.compare(fileName.find_first_not_of(" "),1,".")){ //should be found local
      std::string fullName = wdir+fileName;
      if(checkIfPresent(fullName)){
	cmd.assign("rsync -e ssh -auv "+fullName+" "+_remoteHost+":"+_remoteDir+"/"+fileName);
	return mySystem(cmd);
      }
      else{
	Msg::Fatal("The input file <%s> is not present", fullName.c_str());
	return false;
      }
    }
    else { //should be found remote
      if(!checkIfPresentRemote(fileName)){
	Msg::Fatal("The input file <%s> is not present", fileName.c_str());
	return false;
      }
      else
	return true;
    }
  }
}

bool remoteClient::syncOutputFile(const std::string &wdir, const std::string &fileName){
  std::string cmd;

  if(checkIfPresentRemote(fileName)){
    int pos=fileName.find_first_not_of(" ");
    if(!fileName.compare(pos,1,".")){ // the file must be copied back locally
      cmd.assign("rsync -e ssh -auv "+_remoteHost+":"+_remoteDir+dirSep
		 +fileName.substr(pos,std::string::npos)+" .");
      if(!wdir.empty()) 
	cmd.append(dirSep+wdir);
      return mySystem(cmd);
    }
  }
  else
    return false;
}


// REMOTE INTERFACED Client

bool RemoteInterfacedClient::checkCommandLine(){
  struct stat buf;
  std::string cmd;
  char cbuf [1024];
  FILE *fp;

  if(!isActive()) return true;
  cmd.assign("ssh "+getRemoteHost()+" 'mkdir -p "+getRemoteDir()+"'");
  mySystem(cmd);

  cmd.assign("ssh "+getRemoteHost()+" 'which "+getCommandLine()+"'");
  fp = popen(cmd.c_str(), "r");
  if(fgets(cbuf, 1024, fp) == NULL){
    Msg::Warning("The executable <%s> does not exist", getCommandLine().c_str());
    pclose(fp);
    return false;
  }
  Msg::Info("The executable <%s> exists", getCommandLine().c_str());
  pclose(fp);

  return true;
}

void RemoteInterfacedClient::compute(){
  std::string cmd,rmcmd;
  std::vector<std::string> choices;

  convert();
  Msg::SetOnelabString(getName() + "/Action","compute",false); 

  if(getList("InputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++)
      syncInputFile(getWorkingDir(),choices[i]);
  }

  if(buildRmCommand(rmcmd)){
    cmd.assign("ssh "+getRemoteHost()+" 'cd "+getRemoteDir()+"; "+rmcmd+"'");
    mySystem(cmd);
  }

  cmd.assign("ssh "+getRemoteHost()+" 'cd "+getRemoteDir()+"; "
	     +getCommandLine()+" "+getString("Arguments")+"'");
  mySystem(cmd);

  if(getList("OutputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++)
      syncOutputFile(getWorkingDir(),choices[i]);
  }

  if(getList("PostArray",choices))
    PostArray(choices);
}

// REMOTE ENCAPSULATED Client

std::string RemoteEncapsulatedClient::buildCommandLine(){
  std::string command;
  command.assign("incomp_ssh -f "+getRemoteHost()+" 'cd "+getRemoteDir()+"; "
	         +FixWindowsPath(getCommandLine())+" ");
  if(command.size()){
    std::vector<onelab::string> ps;
    get(ps, getName() + "/Action");
    std::string action = (ps.empty() ? "" : ps[0].getValue());
    get(ps, getName() + "/9CheckCommand");
    std::string checkCommand = (ps.empty() ? "" : ps[0].getValue());
    get(ps, getName() + "/9ComputeCommand");
    std::string computeCommand = (ps.empty() ? "" : ps[0].getValue());

    if(action == "initialize")
      command.append(" " + getSocketSwitch() + " " + getName() + " %s");
    else if(action == "check"){
      command.append(" " + getString("Arguments") + " " + checkCommand + " ");
      command.append(" " + getSocketSwitch() + " " + getName() + " %s");
    }
    else if(action == "compute"){
      command.append(" " + getString("Arguments") + " " + computeCommand + " ");
      command.append(" " + getSocketSwitch() + " " + getName() + " %s");
      //FIXME redirect?
    }
    else
      Msg::Fatal("remoteEncapsulatedClient::run: Unknown Action <%s>", action.c_str());
  }
  return command;
}

bool RemoteEncapsulatedClient::checkCommandLine(){
  struct stat buf;
  std::string cmd;
  char cbuf [1024];
  FILE *fp;

  if(!isActive()) return true;
  cmd.assign("ssh "+getRemoteHost()+" 'mkdir -p "+getRemoteDir()+"'");
  mySystem(cmd);

  cmd.assign("ssh "+getRemoteHost()+" 'which "+getCommandLine()+"'");
  fp = popen(cmd.c_str(), "r");
  if(fgets(cbuf, 1024, fp) == NULL){
    Msg::Warning("The executable <%s> does not exist", getCommandLine().c_str());
    pclose(fp);
    return false;
  }
  Msg::Info("The executable <%s> exists", getCommandLine().c_str());
  pclose(fp);

  Msg::SetOnelabString(getName() + "/Action","initialize",false);
  run();

  return true;
}

void RemoteEncapsulatedClient::analyze(){
  std::string cmd,rmcmd;
  std::vector<std::string> choices;

  Msg::SetOnelabString(getName() + "/Action","check",false);

  if(getList("InputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++)
      syncInputFile(getWorkingDir(),choices[i]);
  }
  run();
}


void RemoteEncapsulatedClient::compute(){
  std::string cmd,rmcmd;
  std::vector<std::string> choices;

  Msg::SetOnelabString(getName() + "/Action","compute",false); 

  if(getList("InputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++)
      syncInputFile(getWorkingDir(),choices[i]);
  }

  if(buildRmCommand(rmcmd)){
    cmd.assign("ssh "+getRemoteHost()+" 'cd "+getRemoteDir()+"; "+rmcmd+"'");
    mySystem(cmd);
  }

  run();

  if(getList("OutputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++)
      syncOutputFile(getWorkingDir(),choices[i]);
  }

  if(getList("PostArray",choices))
    PostArray(choices);
}

// ONELAB additional TOOLS (no access to server in tools)
// options for 'onelab_client'

int getOptions(int argc, char *argv[], std::string &action, std::string &commandLine, std::string &caseName, std::string &clientName, std::string &sockName, int &modelNumber){

  commandLine=argv[0];
  action="compute";
  caseName="untitled";
  int i= 1;
  while(i < argc) {
    if(argv[i][0] == '-') {
      if(!strcmp(argv[i] + 1, "m")) {
	i++;
	modelNumber = (atoi(argv[i]));
        i++;
      }
      else if(!strcmp(argv[i] + 1, "onelab")) {
	i++;
	clientName = argv[i];
        i++;
	sockName = argv[i];
        i++;
      }
      else if(!strcmp(argv[i] + 1, "a")) {
	i++;
	action="check";
      }
      else if(!strcmp(argv[i] + 1, "c")) {
	i++;
	std::cout << "\nONELAB: Present state of the onelab clients\n" << std::endl;
	for(onelab::server::citer it = onelab::server::instance()->firstClient();
	    it != onelab::server::instance()->lastClient(); it++){
	  std::string name= it->second->getName();
	  std::cout << "<" << onelab::server::instance()->getChanged(name) << "> " << name << std::endl;
	}
	action="check";
      }
      else {
	i++;
	printf("Usage: %s [-m num -a -c]\n", argv[0]);
	printf("Options are:\nm      model number\n");
	printf("a      analyze only\n");
	exit(1);
      }
    }
    else{
      caseName=argv[i];
      i++;
    }
  }
}

std::string itoa(const int i){
  std::ostringstream tmp;
  tmp << i ;
  return tmp.str();
}

std::string ftoa(const double x){
  std::stringstream Num;
  Num << x;
  return Num.str();
}

#include <sys/stat.h>		
#include <ctime>
bool checkIfPresent(std::string fileName){
  struct stat buf;
  if (!stat(fileName.c_str(), &buf))
    return true;
  else{
    Msg::Fatal("The file <%s> is not present.",fileName.c_str());
    return false;
  }
}

#include <unistd.h>
#include <sys/types.h>
#if not defined WIN32
#include <pwd.h>
std::string getUserHomedir(){
  struct passwd *pw = getpwuid(getuid());
  std::string str(pw->pw_dir);
  str.append("/");
  return str;
}
#endif

#include <sys/param.h>
std::string getCurrentWorkdir(){
  char path[MAXPATHLEN];
  getcwd(path, MAXPATHLEN);
  std::string str = path;
  return str;
}

std::string sanitize(const std::string &in)
{
  std::string out, forbidden(" ();");
  for(unsigned int i = 0; i < in.size(); i++)
    if ( forbidden.find(in[i]) == std::string::npos)
      out.push_back(in[i]);
  return out;
}
std::string removeBlanks(const std::string &in)
{
  int pos0=in.find_first_not_of(" ");
  int pos=in.find_last_not_of(" ");
  if( (pos0 != std::string::npos) && (pos != std::string::npos))
    return in.substr(pos0,pos-pos0+1);
  else
    return "";
}
bool isPath(const std::string &in)
{
  int pos=in.find_last_not_of(" 0123456789");
  if(in.compare(pos,1,dirSep))
    Msg::Fatal("The argument <%s> is not a valid path (must end with '/')",in.c_str());
  return true;
}

void GmshDisplay(onelab::remoteNetworkClient *loader, std::string fileName, std::vector<std::string> choices){
  if(choices.empty()) return;
#if defined(WIN32)
  std::string cmd = "gmsh.exe ";
#else
  std::string cmd = "gmsh ";
#endif
  cmd.append( fileName + ".geo ");
  for(unsigned int i = 0; i < choices.size(); i++){
    std::string fileName=Msg::GetOnelabString("Arguments/WorkingDir")+choices[i];
    cmd.append(fileName+" ");
    if(Msg::hasGmsh){
      loader->sendMergeFileRequest(fileName);
      Msg::Info("Send merge request <%s>",fileName.c_str());
    }
  }
  if(!Msg::hasGmsh) system(cmd.c_str());
}

std::vector <double> extract_column(const int col, array data){
  std::vector<double> column;
  for ( int i=0; i<data.size(); i++)
    if (  col>0 && col<=data[i].size())
      column.push_back(data[i][col-1]);
    else
      Msg::Fatal("Column number (%d) out of range.",col);
  return column;
}

double find_in_array(int lin, int col, const std::vector <std::vector <double> > &data){
  if ( lin<0 ) {
    lin=data.size();
  }
  if ( lin>=1 && lin<=data.size()){
    if ( col>=1 && col<=data[lin-1].size() )
      return data[lin-1][col-1];
  }
  Msg::Fatal("The value has not been calculated: (%d,%d) out of range",lin,col);
}

array read_array(std::string filename, char sep){
  std::ifstream infile(sanitize(filename).c_str());
  std::vector <std::vector <double> > array;

  int deb,end;
  double temp;
  while (infile){
    std::string s;
    if (!getline( infile, s )) break;
    std::vector <double> record;
    end=0;
    while ( (deb=s.find_first_not_of(" \t\n", end)) != std::string::npos ) {
      end=s.find_first_of(" \t\n",deb);
      temp=atof( s.substr(deb,end).c_str() );
      record.push_back( temp );
    }
    array.push_back( record );
  }
  if (!infile.eof()){
    std::cerr << "Error reading array\n";
  }
  return array;
}

  
