#include "OnelabMessage.h"
#include "OnelabClients.h"
#include "StringUtils.h"
#include <algorithm>

class onelabServer : public GmshServer{
 private:
  localNetworkSolverClient *_client;
 public:
  onelabServer(localNetworkSolverClient *client) : GmshServer(), _client(client) {}
  ~onelabServer() {}
  int SystemCall(const char *str)
  { 
    printf("ONELAB System call(%s)\n", str);
    return system(str); 
  }
  int NonBlockingWait(int socket, double waitint, double timeout)
  { 
    double start = GetTimeInSeconds();
    while(1){
      if(timeout > 0 && GetTimeInSeconds() - start > timeout)
        return 2; // timeout
      if(_client->getPid() < 0)
        return 1; // process has been killed

      // check if there is data (call select with a zero timeout to
      // return immediately, i.e., do polling)
      int ret = Select(0, 0, socket);
      if(ret == 0){ 
        // nothing available: wait at most waitint seconds
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

bool localNetworkSolverClient::run()
{
 new_connection:
  _pid = 0;
  _gmshServer = 0;

  onelabServer *server = new onelabServer(this);
 
#if defined WIN32
  std::string socketName = ":";
#else
  std::string socketName = getUserHomedir() + ".gmshsock";
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
      tmp << GetHostName(); // prepend hostname if only the port number is given
    //tmp << socketName << getId();
    tmp << socketName ;
  }
  sockname = tmp.str();

  std::string command = FixWindowsPath(getCommandLine());
  if(command.size()){
    std::vector<onelab::string> ps;
    get(ps, getName() + "/Action");
    std::string action = (ps.empty() ? "" : ps[0].getValue());
    get(ps, getName() + "/9CheckCommand");
    std::string checkCommand = (ps.empty() ? "" : ps[0].getValue());
    get(ps, getName() + "/9ComputeCommand");
    std::string computeCommand = (ps.empty() ? "" : ps[0].getValue());

    if(action == "initialize")
      command += " ";
    else if(action == "check")
      command += " " + buildArguments() + " " + checkCommand;
    else if(action == "compute")
      command += " " + buildArguments() + " " + computeCommand;
    else
      Msg::Fatal("localNetworkSolverClient::run: Unknown: Unknown Action <%s>", action.c_str());

    // append "-onelab" command line argument
    command += " " + _socketSwitch + " \"" + getName() + "\"";
  }
  else{
    Msg::Info("Listening on socket '%s'", sockname.c_str());
  }

  // std::cout << "FHF sockname=" << sockname.c_str() << std::endl;
  // std::cout << "FHF command=" << command.c_str() << std::endl;

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

  Msg::StatusBar(2, true, "ONELAB: Now running client '%s'...", _name.c_str());
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
        if(type == "number"){
          onelab::number p;
          p.fromChar(message);
          set(p);
        }
        else if(type == "string"){
          onelab::string p;
          p.fromChar(message);
          set(p);
        }
        else
          Msg::Fatal("FIXME query not done for this parameter type: <%s>", message.c_str());
      }
      break;
    case GmshSocket::GMSH_PARAMETER_QUERY:
      {
        std::string version, type, name;
        onelab::parameter::getInfoFromChar(message, version, type, name);
        if(type == "number"){
          std::vector<onelab::number> par;
          get(par, name);
          if(par.size() == 1){
            std::string reply = par[0].toChar();
            server->SendMessage(GmshSocket::GMSH_PARAMETER, reply.size(), &reply[0]);
          }
          else{
            std::string reply = "Parameter (number) " + name + " not found";
            server->SendMessage(GmshSocket::GMSH_INFO, reply.size(), &reply[0]);
          }
        }
        else if(type == "string"){
          std::vector<onelab::string> par;
          get(par, name);
          if(par.size() == 1){
            std::string reply = par[0].toChar();
            server->SendMessage(GmshSocket::GMSH_PARAMETER, reply.size(), &reply[0]);
          }
          else{
            std::string reply = "Parameter (string) " + name + " not found";
            server->SendMessage(GmshSocket::GMSH_INFO, reply.size(), &reply[0]);
          }
        }
        else
          Msg::Fatal("FIXME query not done for this parameter type: <%s>", message.c_str());
      }
      break;
    case GmshSocket::GMSH_PARAM_QUERY_ALL:
      {
        std::string version, type, name, reply;
        onelab::parameter::getInfoFromChar(message, version, type, name);
	if(type == "number"){
	  std::vector<onelab::number> numbers;
	  get(numbers, "");
	  for(std::vector<onelab::number>::iterator it = numbers.begin(); it != numbers.end(); it++){
	    reply = (*it).toChar();
	    server->SendMessage(GmshSocket::GMSH_PARAM_QUERY_ALL, reply.size(), &reply[0]);
	  }
	  server->SendMessage(GmshSocket::GMSH_PARAM_QUERY_END, 0, NULL);
	}
	else if(type == "string"){
	  std::vector<onelab::string> strings;
	  get(strings, "");
	  for(std::vector<onelab::string>::iterator it = strings.begin(); it != strings.end(); it++){
	    reply = (*it).toChar();
	    server->SendMessage(GmshSocket::GMSH_PARAM_QUERY_ALL, reply.size(), &reply[0]);
	  }
	  server->SendMessage(GmshSocket::GMSH_PARAM_QUERY_END, 0, NULL);
	}
        else
          Msg::Fatal("FIXME query not done for this parameter type: <%s>", message.c_str());
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
      SystemCall("gmsh "+ message+" &");
      break;
    default:
      Msg::Warning("Received unknown message type (%d)", type);
      break;
    }
  }

  server->Shutdown();
  delete server;
  Msg::StatusBar(2, true, "ONELAB: Done running '%s'", _name.c_str());
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

// PROMPTUSER

int PromptUser::getVerbosity(){
  std::vector<onelab::number> numbers;
  get(numbers,"VERBOSITY");
  if (numbers.size())
    return numbers[0].getValue();
  else
    return 0;
}
void PromptUser::setVerbosity(const int ival){
  onelab::number number;
  number.setName("VERBOSITY");
  number.setValue(ival);
  set(number);
}

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
  std::string db = "ONELAB parameter space: size=" + itoa(onelab::server::instance()->getNumParameters()) + "\n";
  db.append(onelab::server::instance()->toChar());
  for(unsigned int i = 0; i < db.size(); i++)
    if(db[i] == onelab::parameter::charSep()) db[i] = '|';
  return db.c_str();
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

bool PromptUser::menu(std::string commandLine, std::string fileName, int modelNumber) { 
  int choice, counter1=0, counter2=0;
  std::string answ;
  std::string name;
  onelab::number x;
  onelab::string y;
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;

  std::string clientName="loadedMetaModel";
  EncapsulatedClient *loadedSolver = new  EncapsulatedClient(clientName,commandLine);
  //setString("Arguments/FileName",fileName);
  addStringChoice(clientName + "/InputFiles",fileName);

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
      // std::cout << "\nONELAB: Present state of the parameter space\n" << std::endl;
      // std::string db = onelab::server::instance()->toChar();
      // for(unsigned int i = 0; i < db.size(); i++)
      // 	if(db[i] == onelab::parameter::charSep()) db[i] = '|';
      // std::cout << db.c_str();
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

// LOCALSOLVERCLIENT

bool localSolverClient::controlPath(){
  std::string commandLine = Msg::GetOnelabString(getName()+"/Path");
  if(commandLine.size())
    setCommandLine(commandLine);
  if(getCommandLine().empty()){
    if(Msg::hasGmsh) {// exits metamodel and restores control to the onelab window
      Msg::Error("The path to client <%s> is undefined.", getName().c_str());
      std::cout << "\n\nEnter the path to <" << getName() << "> in the ONELAB window.\n\n" << std::endl;
      //Msg::SetOnelabAttributeString(getName()+"/Path","Highlight","true");
      return false;
    }
    else{ // asks the user in console mode
      std::cout << "\nONELAB:Enter the path on your system to the executable file of <" << getName() << ">" << std::endl;
      std::string path;
      std::getline (std::cin,path);
      setCommandLine(path);
      return true;
    }
  }
  else{ //initialize
    Msg::SetOnelabString(getName() + "/Action","initialize",false);
    run();
  }
  return true;
}

const std::string localSolverClient::getLineOptions(){
  std::vector<onelab::string> strings;
  std::string paramName(getName()+"/LineOptions");
  get(strings,paramName);
  if(strings.size())
    return strings[0].getValue();
  else
    return "";
}
const std::string localSolverClient::getPreLineOptions(){
  std::vector<onelab::string> strings;
  std::string paramName(getName()+"/PreLineOptions");
  get(strings,paramName);
  if(strings.size())
    return strings[0].getValue();
  else
    return "";
}

const std::vector<std::string> localSolverClient::getInputFiles(){
  std::vector<onelab::string> strings, choices;
  std::string fileName;
  std::string paramName(getName()+"/InputFiles");
  get(strings,paramName);
  if(strings.size())
    return strings[0].getChoices();
  else
    Msg::Fatal("ONELAB: parameter <%s> undefined",paramName.c_str());
}

const std::string localSolverClient::buildArguments(){
  int pos;
  std::vector<onelab::string> strings;
  std::string args,filename;

  std::vector<std::string> choices = getInputFiles();
  args.assign(getPreLineOptions()+" ");
  for(unsigned int i = 0; i < choices.size(); i++)
      args.append(choices[i].substr(0,choices[i].find(olkey::extension))+" ");
  args.append(getLineOptions());
  //std::cout << "args=<" << args << ">" << std::endl;
  return args;
}

// METAMODEL

void MetaModel::savePathes(const std::string fileName){ // save client pathes
  std::string fileNameSave = fileName + olkey::extension + ".save";
  std::ofstream outfile(fileNameSave.c_str()); 
  if (outfile.is_open())
    for(citer it = _clients.begin(); it != _clients.end(); it++)
      outfile << (*it)->toChar();
  else
    Msg::Fatal("The file <%s> cannot be opened",fileNameSave.c_str());
  outfile.close();
}

bool MetaModel::checkPathes(){
  bool allDefined=true;
  for(citer it = _clients.begin(); it != _clients.end(); it++){
    allDefined = allDefined && (*it)->controlPath();
  }
  savePathes(genericNameFromArgs);
  return allDefined;
}

void MetaModel::initialize()
{
  Msg::Info("Metamodel::initialize <%s>",getName().c_str());
  Msg::SetOnelabString(clientName + "/9CheckCommand","-a",false);
  Msg::SetOnelabNumber(clientName + "/UseCommandLine",1,false);
  Msg::SetOnelabNumber(clientName + "/Initialized",1,false);
}

std::string MetaModel::resolveGetVal(std::string line) {
  //std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  std::vector<std::string> arguments;
  std::string buff;
  int pos,pos0,cursor;

  cursor=0;
  while ( (pos=line.find(olkey::getValue,cursor)) != std::string::npos){
    pos0=pos; // for further use
    cursor = pos+olkey::getValue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos+1-cursor),arguments)<1)
      Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::getValue.c_str(),line.c_str());
    get(strings,arguments[0]);
    if (strings.size())
      buff.assign(strings[0].getValue());
    else
      Msg::Fatal("ONELAB unknown variable: %s",arguments[0].c_str());
    line.replace(pos0,pos-pos0,buff); 
    cursor=pos0+buff.length();
  }
  return line;
}

void MetaModel::registerClient(const std::string name, const std::string type, 
			       const std::string path) {
  localSolverClient *c;
  Msg::Info("ONELAB: initialize client <%s>", name.c_str());
  if(!type.compare("encapsulated"))  
    c = new EncapsulatedClient(name,path);
  else if(!type.compare("interfaced"))
    c = new InterfacedClient(name,path);
  else
    Msg::Fatal("ONELAB: unknown client type <%s>",type.c_str());
  _clients.push_back(c); 
  if(path.size()){
    //cannot initialize a client without a path
    Msg::SetOnelabString(name + "/Action","initialize",false);
    c->run();
  }
}

void MetaModel::simpleCheck()
{
  for(citer it = _clients.begin(); it != _clients.end(); it++){
    if((*it)->getActive()){
	Msg::SetOnelabString((*it)->getName() + "/Action","check",false);
	(*it)->analyze();
      }
  }
}

void MetaModel::simpleCompute()
{
  for(citer it = _clients.begin(); it != _clients.end(); it++){
    if((*it)->getActive()){
	Msg::SetOnelabString((*it)->getName() + "/Action","compute",false);
	(*it)->compute();
      }
  }
}

void MetaModel::PostArray(std::vector<std::string> choices)
{
  int nb=0;
  onelab::number o;
  while( 4*(nb+1) <= choices.size()){
    //std::cout << "Nb Choices" << choices.size() << std::endl;
    int lin= atof(choices[4*nb+1].c_str());
    int col= atof(choices[4*nb+2].c_str());
    double val=find_in_array(lin,col,read_array(choices[4*nb],' '));
    // o.setName(choices[4*nb+3]);
    // o.setValue(val);
    // set(o);
    Msg::AddOnelabNumberChoice(choices[4*nb+3],val);
    Msg::Info("PostArray <%s>=%e",choices[4*nb+3].c_str(),val);
    nb++;
  }
}

// INTERFACED client

std::string  InterfacedClient::toChar() {
  std::ostringstream sstream;
  if(getCommandLine().size()){
    sstream << olkey::client << " " 
	    << getName() << "." << "Path("
	    << getCommandLine() << ");\n";
  }
  return sstream.str();
}

void InterfacedClient::analyze() { 
  std::vector<onelab::string> strings;
  Msg::SetOnelabString(getName() + "/Action","check",false);// a titre indicatif

  get(strings,getName()+"/InputFiles");
  if(strings.size()){
    std::vector<std::string> choices=strings[0].getChoices();
    for(unsigned int i = 0; i < choices.size(); i++){
      std::string ifilename = choices[i];
      if(ifilename.find(olkey::extension)!=std::string::npos){
	checkIfPresent(ifilename);
	parse_onefile(ifilename); // recursive
      }
    }
  }
}

void InterfacedClient::convert() {
  int pos;
  std::vector<std::string> choices;
  if(Msg::GetOnelabChoices(getName()+"/InputFiles",choices)){
    for(unsigned int i = 0; i < choices.size(); i++){
      std::string ifilename = choices[i];
      checkIfPresent(ifilename);
      if((pos=ifilename.find(olkey::extension))!=std::string::npos){
	std::string ofilename = ifilename.substr(0,pos);  // remove extension
	std::ofstream outfile(ofilename.c_str());
	if (outfile.is_open())
	  convert_onefile(ifilename,outfile);
	else
	  Msg::Fatal("The file <%s> cannot be opened",ofilename.c_str());
	outfile.close();
	//checkIfModified(ofilename); not really useful if opening => touching
      }
    }
  }
}

void InterfacedClient::compute() { 
  convert();
  Msg::SetOnelabString(getName() + "/Action","compute",false); // a titre indicatif

  std::string commandLine = getCommandLine() + " " ;
  commandLine.append(buildArguments());
  //commandLine.append(" &> " + _name + ".log");
  Msg::Info("Client %s launched",_name.c_str());
  if ( int error = systemCall(commandLine.c_str()))
    Msg::Error("Client %s returned error %d",_name.c_str(),error);
  std::vector<std::string> choices;
  if(Msg::GetOnelabChoices(getName()+"/OutputFiles",choices))
    checkIfModified(choices);
  Msg::Info("Client %s completed",_name.c_str());
}

// ENCAPSULATED Client

std::string EncapsulatedClient::toChar(){
  std::ostringstream sstream;
  if(getCommandLine().size()){
    sstream << olkey::client << " " 
	    << getName() << "." << "Path(" 
	    << getCommandLine() << ");\n";
  }
  return sstream.str();
}

void EncapsulatedClient::analyze() {
  set(onelab::string(getName()+"/Action", "check"));
  run();
}

void EncapsulatedClient::compute() {
  set(onelab::string(getName()+"/Action", "compute"));
  run();
}


// ONELAB additional TOOLS (no access to server in tools)

// utilisé pour le main() des métamodèles
int getOptions(int argc, char *argv[], std::string &action, std::string &commandLine, std::string &fileName, std::string &clientName, std::string &sockName, int &modelNumber){

  commandLine=argv[0];
  action="compute";
  fileName="untitled";
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
      fileName=argv[i];
      i++;
    }
  }
}

static std::string getNextToken(const std::string &msg,
				std::string::size_type &first){
  if(first == std::string::npos) return "";
  std::string::size_type last = msg.find_first_of(charSep(), first);
  std::string next = msg.substr(first, last - first);
  first = (last == std::string::npos) ? last : last + 1;
  return next;
}
 
std::string itoa(const int i){
  std::ostringstream tmp;
  tmp << i ;
  return tmp.str();
}

int onelab_step;
int newStep(){
  if (onelab_step==0){
    system("echo 'start' > onelab.log");
    system("touch onelab.start; touch onelab.progress");
  }
  system("find . -newer onelab.progress > onelab.modified");
  system("touch onelab.progress");
  onelab_step++;
}
int getStep(){
  return onelab_step;
}

std::string sanitize(const std::string &in)
{
  std::string out, forbidden(" ();");
  for(unsigned int i = 0; i < in.size(); i++)
    if ( forbidden.find(in[i]) == std::string::npos)
      out.push_back(in[i]);
  return out;
}
int enclosed(const std::string &in, std::vector<std::string> &arguments){
  int pos, cursor;
  arguments.resize(0);
  cursor=0;
  if ( (pos=in.find("(",cursor)) == std::string::npos )
     Msg::Fatal("ONELAB: syntax error: <%s>",in.c_str());

  unsigned int count=1;
  pos++; // skips '('
  cursor = pos; 
  do{
    if(in[pos]=='(') count++;
    if(in[pos]==')') count--;
    if(in[pos]==',') {
      arguments.push_back(in.substr(cursor,pos-cursor));
      if(count!=1)
	Msg::Fatal("ONELAB: syntax error: <%s>",in.c_str());
      cursor=pos+1; // skips ','
    }
    pos++;
  } while( count && (pos!=std::string::npos) ); // find closing brace
  if(count)
     Msg::Fatal("ONELAB: syntax error: <%s>",in.c_str());
  else
    arguments.push_back(in.substr(cursor,pos-1-cursor));
  return arguments.size();
}
int extract(const std::string &in, std::string &paramName, std::string &action, std::vector<std::string> &arguments){
  // syntax: paramName.action( arg1, arg2, ... )
  int pos, cursor,NumArg=0;
  cursor=0;
  if ( (pos=in.find(".",cursor)) == std::string::npos )
     Msg::Fatal("ONELAB: syntax error: <%s>",in.c_str());
  else
    paramName.assign(sanitize(in.substr(cursor,pos-cursor)));
  cursor = pos+1; // skips '.'
  if ( (pos=in.find("(",cursor)) == std::string::npos )
     Msg::Fatal("ONELAB: syntax error: <%s>",in.c_str());
  else
    action.assign(sanitize(in.substr(cursor,pos-cursor)));
  cursor = pos;
  unsigned int count=0;
  do{
    if(in[pos]=='(') count++;
    if(in[pos]==')') count--;
    pos++;
  } while(count && (pos!=std::string::npos) ); // find closing brace
  if(count)
     Msg::Fatal("ONELAB: syntax error: %s",in.c_str());
  else
    NumArg = enclosed(in.substr(cursor,pos-cursor),arguments);
  // std::cout << "paramName=<" << paramName << ">" << std::endl;
  // std::cout << "arguments=<" << in.substr(cursor,pos+1-cursor) << ">" << std::endl;
  return NumArg;
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

#include <sys/stat.h>		
#include <ctime>

bool fileExist(std::string filename){
  struct stat buf;
  if(!stat(filename.c_str(), &buf)){
    std::string cmd = "touch " + filename;
    system(cmd.c_str());
    return true;
  }
  else
    return false;
}

bool checkIfPresent(std::string filename){
  struct stat buf;
  if (!stat(filename.c_str(), &buf))
    return true;
  else{
    Msg::Fatal("The file %s is not present",filename.c_str());
    return false;
  }
}

bool checkIfModified(std::vector<std::string> filenames){
  for(unsigned i=0; i<filenames.size(); i++)
    checkIfModified(filenames[i]);
  return true;
}
bool checkIfModified(std::string filename){
  struct stat buf1,buf2;
  if (stat("onelab.start", &buf1))
    Msg::Fatal("The file %s does not exist.","onelab.start");
  if (stat(filename.c_str(), &buf2))
    Msg::Fatal("The file %s does not exist.",filename.c_str());
  if (difftime(buf1.st_mtime, buf2.st_mtime) > 0)
    Msg::Fatal("The file %s has not been modified.",filename.c_str());
  return true;
}

int systemCall(std::string cmd){
  printf("ONELAB: System call <%s>\n", cmd.c_str());
  int err=system(cmd.c_str());
  printf("ONELAB: System call <%s> returns <%d>\n", cmd.c_str(),err);
  return err;
}

void GmshDisplay(onelab::remoteNetworkClient *loader, std::string fileName, std::vector<std::string> choices){
  if(choices.empty()) return;
  std::string cmd = "gmsh " + fileName + ".geo ";
  for(unsigned int i = 0; i < choices.size(); i++){
    cmd.append(choices[i]+" ");
    checkIfModified(choices[i]);
    if(Msg::hasGmsh){
      loader->sendMergeFileRequest(choices[i]);
      Msg::Info("Send merge request <%s>",choices[i].c_str());
    }
  }
  if(!Msg::hasGmsh) systemCall(cmd);
}
void GmshDisplay(onelab::remoteNetworkClient *loader, std::string modelName, std::string fileName){
  checkIfModified(fileName);
  if(loader)
    loader->sendMergeFileRequest(fileName);
  else {
    std::string cmd= "gmsh " + modelName + ".geo " + fileName;
    systemCall(cmd);
  }
}

void appendOption(std::string &str, const std::string &what, const int val){
  if(val){ // assumes val=0 is the default
    std::stringstream Num;
    Num << val;
    str.append( what + " " + Num.str() + " ");
  }
}
void appendOption(std::string &str, const std::string &what){
  str.append( what + " ");
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

double find_in_array(const int lin, const int col, const std::vector <std::vector <double> > &data){
  if ( lin>=1 && lin<=data.size())
    if (  col>=1 && col<=data[lin-1].size())
      return data[lin-1][col-1];
  Msg::Fatal("The value has not been calculated: (%d,%d) out of range",lin,col);
}

array read_array(std::string filename, char sep){
  std::ifstream infile(sanitize(filename).c_str());
  std::vector <std::vector <double> > array;

  while (infile){
    std::string s;
    if (!getline( infile, s )) break;
    std::istringstream ss( s );
    std::vector <double> record;
    while (ss){
      std::string s;
      if (!getline( ss, s, sep )) break;
      if ( s.size() ){
	//std::cout << "Read=<" << s << ">" << std::endl;
	record.push_back( atof( s.c_str() ));
      }
    }
    array.push_back( record );
  }
  if (!infile.eof()){
    std::cerr << "Error reading array\n";
  }
  return array;
}

  
