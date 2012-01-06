#include "OnelabClients.h"
#include "OnelabMessage.h"

class onelabServer : public GmshServer{
 private:
  onelab::localNetworkClient *_client;
 public:
  onelabServer(onelab::localNetworkClient *client) : GmshServer(), _client(client) {}
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

bool onelab::localNetworkClient::run(const std::string &what)
{
 new_connection:
  _pid = 0;
  _gmshServer = 0;

  onelabServer *server = new onelabServer(this);
 
  //std::string socketName = ":";
  std::string socketName = getUserHomedir() + ".gmshsock";
  std::string sockname;
  std::ostringstream tmp;
  if(!strstr(socketName.c_str(), ":")){
    // Unix socket
    tmp << socketName << getId();
    //sockname = FixWindowsPath(tmp.str());
  }
  else{
    // TCP/IP socket
    if(socketName.size() && socketName[0] == ':')
      tmp << GetHostName(); // prepend hostname if only the port number is given
    //tmp << socketName << getId();
    tmp << socketName ;
  }
  sockname = tmp.str();

  //std::string command = FixWindowsPath(_commandLine);
  std::string command = _commandLine;
  if(command.size()){
    command += " " + what + " " + _socketSwitch + " ";
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
          Msg::Error("FIXME not done for this parameter type");
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
          Msg::Fatal("FIXME query not done for this parameter type");
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
          Msg::Fatal("FIXME query not done for this parameter type: %s", message.c_str());
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

bool onelab::localNetworkClient::kill()
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
// int PromptUser::getInteractivity(){
//   std::vector<onelab::number> numbers;
//   get(numbers,"INTERACTIVITY");
//   if (numbers.size())
//     return numbers[0].getValue();
//   else
//     return 1;
// }
// void PromptUser::setInteractivity(const int ival){
//   onelab::number number;
//   number.setName("INTERACTIVITY");
//   number.setValue(ival);
//   set(number);
// }
void PromptUser::setNumber(const std::string paramName, const double val, const std::string &str){
  onelab::number number;
  number.setName(paramName);
  number.setValue(val);
  number.setHelp(str);
  set(number);
}
double PromptUser::getNumber(const std::string paramName){
  std::vector<onelab::number> numbers;
  get(numbers,paramName);
  if (numbers.size())
    return numbers[0].getValue();
  else
    Msg::Fatal("Unknown parameter %s",paramName.c_str());
}
void PromptUser::setString(const std::string paramName, const std::string &val, const std::string &str){
  onelab::string string;
  string.setName(paramName);
  string.setValue(val);
  string.setHelp(str);
  set(string);
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
  std::string db = "ONELAB: parameter space state:\n";
  db.append(onelab::server::instance()->toChar());
  for(unsigned int i = 0; i < db.size(); i++)
    if(db[i] == onelab::parameter::charSep()) db[i] = '|';
  return db.c_str();
}

bool PromptUser::menu(std::string options, std::string modelName) { 
  int choice, counter1=0,counter2=0;
  std::string answ;
  std::string name;
  onelab::number x;
  onelab::string y;
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;

  MetaModel *metaModel = new  MetaModel(modelName);

  do {
    std::cout << "\nONELAB: menu" << std::endl ;
    std::cout << " 1- View parameter space\n 2- Set a value\n 3- List modified files\n 4- Analyze\n 5- Compute\n 6- Quit metamodel" << std::endl;

    choice=0;
    std::string mystr;
    while( (choice<1 || choice>6) && ++counter1<10 ) {
      std::cout << "\nONELAB: your choice? "; 
      std::getline (std::cin, mystr);
      std::stringstream myStream(mystr);
      if (myStream >> choice) break;
      std::cout << "Invalid choice" << std::endl;
    }
    std::cout << "Your choice is <" << choice << ">" << std::endl;

    if (choice==1){
      std::cout << "\nONELAB: Present state of the parameter space\n" << std::endl;
      std::string db = onelab::server::instance()->toChar();
      for(unsigned int i = 0; i < db.size(); i++)
	if(db[i] == onelab::parameter::charSep()) db[i] = '|';
      std::cout << db.c_str();
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
    else if (choice==4){
      metaModel->analyze(options+" -a","");
      choice=0;
    }
    else if (choice==5){
      metaModel->run(options,"");
      choice=0;
    }
    else if (choice==6)
      exit(1);
    else
      choice=0;
  } while(!choice && ++counter2<20);
}



std::string sanitize(const std::string &in)
{
  std::string out, forbidden(" ();");
  for(unsigned int i = 0; i < in.size(); i++)
    if ( forbidden.find(in[i]) == std::string::npos)
      out.push_back(in[i]);
  return out;
}
int extract(const std::string &in, std::string &paramName, std::string &action, std::string &value){
  int pos, cursor;
  cursor=0;
  if ( (pos=in.find("."))  == std::string::npos )
     Msg::Fatal("Onelab syntax error: %s",in.c_str());
  else
    paramName.assign(sanitize(in.substr(cursor,pos-cursor)));
  cursor = pos+1;
  if ( (pos=in.find("("),cursor) == std::string::npos )
     Msg::Fatal("Onelab syntax error: %s",in.c_str());
  else
    action.assign(sanitize(in.substr(cursor,pos-cursor)));
  cursor = pos+1;
  if ( (pos=in.find(")"),cursor) == std::string::npos )
     Msg::Fatal("Onelab syntax error: %s",in.c_str());
  else
    value.assign(in.substr(cursor,pos-cursor));
  //std::cout << "param=<" << paramName << "> action=<" << action << "> value=<" << value << ">"<< std::endl;
}
int enclosed(const std::string &in, std::string &out){
  int pos, cursor;
  cursor=0;
  if ( (pos=in.find("("),cursor) == std::string::npos )
     Msg::Fatal("Onelab syntax error: %s",in.c_str());
  cursor = pos+1;
  if ( (pos=in.find(")"),cursor) == std::string::npos )
     Msg::Fatal("Onelab syntax error: %s",in.c_str());
  else
    out.assign(sanitize(in.substr(cursor,pos-cursor)));
  //std::cout << "enclosed=<" << out << ">"<< std::endl;
}

bool InterfacedClient::analyze_oneline(std::string line, std::ifstream &infile) { 
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  int pos0,pos,cursor;
  char sep=';';
  std::string buff;
  std::string onelab("onelab"), number("onelab.number"), include("onelab.include"), iftrue("onelab.iftrue");

  if ( (pos=line.find(number)) != std::string::npos) {// onelab.number
    cursor = pos+number.length();
    while ( (pos=line.find_first_of(sep,cursor)) != std::string::npos){
      //std::cout << line.substr(cursor,pos-cursor) << std::endl;
      std::string name, action, value;
      extract(line.substr(cursor,pos-cursor),name,action,value);
      std::vector<onelab::number> ps;
      get(ps, name);

      if(!ps.size()){ // param does not exist
	if(!action.compare("Create") || !action.compare("Default")){ 
	  ps.resize(1);
	  ps[0].setName(name);
	  ps[0].setValue(atof(value.c_str()));
	}
	else{
	  Msg::Fatal("ONELAB: parameter '%s' has not been created",name.c_str());
	} 
      }
      else{ // param already created
	if (!action.compare("Min"))
	  ps[0].setMin(atof(value.c_str()));
	else if (!action.compare("Max"))
	  ps[0].setMax(atof(value.c_str()));
	else if (!action.compare("Help"))
	  ps[0].setHelp(value.c_str());
	else if(!action.compare("Create") || !action.compare("Default")) {}
	else
	  Msg::Fatal("ONELAB unknown action: %s",action.c_str());
      }
      set(ps[0]); 
      //std::cout << ps[0]->toChar() << std::endl;
      cursor=pos+1;
    }
  }
  else if ( (pos=line.find(iftrue)) != std::string::npos) {// onelab.iftrue
    cursor = pos+iftrue.length();
    pos=line.find_first_of(')',cursor)+1;
    std::string boolParam;
    enclosed(line.substr(cursor,pos-cursor),boolParam);
    std::cout << "iftrue " << boolParam << std::endl;
    
    get(numbers,boolParam);
    if (numbers.size()){
      bool condition = (bool) numbers[0].getValue();
      if (!analyze_ifstatement(infile,condition))
	Msg::Fatal("ONELAB misformed onelab.iftrue statement: %s",boolParam.c_str());
    }
    else
      Msg::Fatal("ONELAB unknown boolean parameter: <%s>",boolParam.c_str());
  }
  else if ( (pos=line.find(include)) != std::string::npos) {// onelab.include
    cursor = pos+include.length();
    pos=line.find_first_of(')',cursor)+1;
    std::string fileName;
    enclosed(line.substr(cursor,pos-cursor),fileName);
    analyze_onefile(fileName);
  }
}

bool InterfacedClient::analyze_onefile(std::string ifilename) { 
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  int pos0,pos,cursor;
  char sep=';';
  std::string line,buff;
  std::string onelab("onelab"), number("onelab.number"), include("onelab.include");
  std::ifstream infile(ifilename.c_str());

  if (infile.is_open()){
    while ( infile.good() ) {
      getline (infile,line);
      analyze_oneline(line,infile);
    }
    infile.close();
  }
  else
    Msg::Fatal("The file %s cannot be opened",ifilename.c_str());
} 

bool InterfacedClient::analyze_ifstatement(std::ifstream &infile, bool condition) { 
  int pos;
  std::string line;
  std::string iftrue("onelab.iftrue"), olelse("onelab.else"), olendif("onelab.endif");

  bool trueclause=true, loopend=false;
  while ( infile.good() && !loopend) {
    getline (infile,line);
    if ( (pos=line.find(olelse)) != std::string::npos) 
      trueclause=false;
    else if ( (pos=line.find(olendif)) != std::string::npos) 
      loopend=true;
    else if ( !(trueclause ^ condition) ) // xor bitwise operator
      analyze_oneline(line,infile);
  }
  return loopend;
} 

bool InterfacedClient::convert_ifstatement(std::ifstream &infile, std::ofstream &outfile, bool condition) { 
  int pos;
  std::string line;
  std::string iftrue("onelab.iftrue"), olelse("onelab.else"), olendif("onelab.endif");

  bool trueclause=true, loopend=false;
  while ( infile.good() && !loopend) {
    getline (infile,line);
    if ( (pos=line.find(olelse)) != std::string::npos) 
      trueclause=false;
    else if ( (pos=line.find(olendif)) != std::string::npos) 
      loopend=true;
    else if ( !(trueclause ^ condition) ) // xor bitwise operator
      convert_oneline(line,infile,outfile);
  }
  return loopend;
} 

bool InterfacedClient::convert_oneline(std::string line, std::ifstream &infile, std::ofstream &outfile) { 
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  int pos0,pos,cursor;
  char sep=';';
  std::string buff;
  std::string onelab("onelab"), number("onelab.number"), include("onelab.include"), getValue("onelab.getValue");
  std::string iftrue("onelab.iftrue");

  if ( (pos=line.find(onelab)) == std::string::npos) // not a onelab line
    outfile << line << std::endl; 
  else{ 
    if ( (pos=line.find(number)) != std::string::npos) {// onelab.number
      //skip the line
    }
    else if ( (pos=line.find(include)) != std::string::npos) {// onelab.include
      cursor = pos+include.length();
      pos=line.find_first_of(')',cursor)+1;
      std::string fileName;
      enclosed(line.substr(cursor,pos-cursor),fileName);
      convert_onefile(fileName,outfile);
    }
    else if ( (pos=line.find(iftrue)) != std::string::npos) {// onelab.iftrue
      cursor = pos+iftrue.length();
      pos=line.find_first_of(')',cursor)+1;
      std::string boolParam;
      enclosed(line.substr(cursor,pos-cursor),boolParam);
      get(numbers,boolParam);
      if (numbers.size()){
	bool condition = (bool) numbers[0].getValue();
	if (!convert_ifstatement(infile,outfile,condition))
	  Msg::Fatal("ONELAB misformed onelab.iftrue statement: %s",boolParam.c_str());
      }
      else
	Msg::Fatal("ONELAB unknown boolean parameter: %s",boolParam.c_str());
    }
    else if ( (pos=line.find(getValue)) != std::string::npos) {
      // onelab.getValue, possibly several times in the line
      cursor=0;
      while ( (pos=line.find(getValue,cursor)) != std::string::npos){
	pos0=pos; // for further use
	cursor = pos+getValue.length();
	pos=line.find_first_of(')',cursor)+1;
	std::string paramName;
	enclosed(line.substr(cursor,pos-cursor),paramName);
	//std::cout << line.substr(pos0,pos) << std::endl;

	get(numbers,paramName);
	if (numbers.size()){
	  std::stringstream Num;
	  Num << numbers[0].getValue();
	  buff.assign(Num.str());
	}
	else{
	  get(strings,paramName);
	  if (strings.size())
	    buff.assign(strings[0].getValue());
	  else
	    Msg::Fatal("ONELAB unknown variable: %s",paramName.c_str());
	}
	line.replace(pos0,pos-pos0,buff); 
	cursor=pos0+buff.length();
      }
      outfile << line << std::endl; 
    }
    else
      Msg::Fatal("ONELAB syntax error: %s",line.c_str());
  }
}

bool InterfacedClient::convert_onefile(std::string ifilename, std::ofstream &outfile) { 
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  int pos0,pos,cursor;
  char sep=';';
  std::string line,buff;
  std::string onelab("onelab"), number("onelab.number"), include("onelab.include"), getValue("onelab.getValue");
  std::ifstream infile(ifilename.c_str());

  if (infile.is_open()){
    while ( infile.good() ) {
      getline (infile,line);
      convert_oneline(line,infile,outfile);
    }
    infile.close();
  }
  else
    Msg::Fatal("The file %s cannot be opened",ifilename.c_str());
}

bool InterfacedClient::analyze(const std::string options, const std::string modelName) { 
  std::string ifilename = modelName+ _extension + "_onelab";
  analyze_onefile(ifilename); // recursive
}

bool InterfacedClient::convert(std::string modelName) { 
  std::string ifilename = modelName+ _extension + "_onelab";
  std::string ofilename = modelName+ _extension ;
  std::ofstream outfile(ofilename.c_str());

  if (outfile.is_open())
    convert_onefile(ifilename,outfile);
  else
    Msg::Fatal("The file %s cannot be opened",ofilename.c_str());
  outfile.close();
}

bool InterfacedClient::run(const std::string options, const std::string modelName) { 
  convert(modelName);
  std::string commandLine = _commandLine + " " + modelName + _extension;
  commandLine.append(" " +options);
  //commandLine.append(" &> " + _name + ".log");
  Msg::Info("Client %s launched",_name.c_str());
  std::cout << "Commandline:" << commandLine.c_str() << std::endl;
  if ( int error = system(commandLine.c_str())) { 
    Msg::Error("Client %s returned error %d",_name.c_str(),error);
  }
  Msg::Info("Client %s completed",_name.c_str());
  return true;
}

bool EncapsulatedTest::analyze(const std::string options, const std::string modelName) {
  // onelab::server::citer it= onelab::server::instance()->findClient(getName());
  // onelab::client *c = it->second;
  // c->run("");
  run(options,modelName);
  return true;
}

bool EncapsulatedTest::run(const std::string options, const std::string modelName) {
  // onelab::server::citer it= onelab::server::instance()->findClient(getName());
  // onelab::client *c = it->second;
  // std::string commandLine = modelName + options;
  // c->run("");
  run(options,modelName);
  return true;
}

int EncapsulatedGmsh::getVerbosity(){
  std::vector<onelab::number> numbers;
  get(numbers,"VERBOSITY");
  if (numbers.size())
    return numbers[0].getValue();
  else
    return 0;
}
bool EncapsulatedGmsh::analyze(const std::string options, const std::string modelName) {
  onelab::server::citer it= onelab::server::instance()->findClient(getName());
  onelab::client *c = it->second;
  std::string commandLine = modelName + ".geo - " + options;
  appendOption(commandLine,"-v",getVerbosity());
  c->run(commandLine);
  return true;
}

bool EncapsulatedGmsh::run(const std::string options, const std::string modelName) {
  onelab::server::citer it= onelab::server::instance()->findClient(getName());
  onelab::client *c = it->second;
  std::string commandLine = modelName + ".geo " + options ;
  appendOption(commandLine,"-v",getVerbosity());
  c->run(commandLine);
  return true;
}

// int EncapsulatedGetdp::appendVerbosity(std::string &str){
//   std::vector<onelab::number> numbers;
//   int verb=0;
//   get(numbers,"VERBOSITY");
//   if (numbers.size())
//     verb = numbers[0].getValue();
//   std::stringstream Num;
//   Num << verb;
//   str.append(" -v " + Num.str() + " ");
//   return verb;
// }
// int EncapsulatedGetdp::appendResolution(std::string &str){
//   std::vector<onelab::string> strings;
//   get(strings,"GetDP/1resolution");
//   if (strings.size())
//     str.append(" -sol " + strings[0].getValue() + " ");
//   else
//     Msg::Fatal("Resolution <GetDP/1resolution> not defined");
//   return strings.size();
// }
// int EncapsulatedGetdp::appendPostpro(std::string &str){
//   std::vector<onelab::string> strings;
//   get(strings,"GetDP/2Post-Operation");
//   if (strings.size())
//     str.append(" -pos " + strings[0].getValue() + " ");
//   else
//     Msg::Fatal("Resolution <GetDP/2Post-Operation> not defined");
//   return strings.size();
// }
int EncapsulatedGetdp::getVerbosity(){
  std::vector<onelab::number> numbers;
  get(numbers,"VERBOSITY");
  if (numbers.size())
    return numbers[0].getValue();
  else
    return 0;
}
bool EncapsulatedGetdp::analyze(const std::string options, const std::string modelName) {
  onelab::server::citer it= onelab::server::instance()->findClient(getName());
  onelab::client *c = it->second;
  std::string commandLine = modelName + " " + options  ;
  appendOption(commandLine,"-v",getVerbosity());
  c->run(modelName);
  return 1;
}
bool EncapsulatedGetdp::run(const std::string options, const std::string modelName) {
  onelab::server::citer it= onelab::server::instance()->findClient(getName());
  onelab::client *c = it->second;
  std::string commandLine = modelName + " " + options;
  appendOption(commandLine,"-v",getVerbosity());
  c->run(commandLine);
  return true;
}
// bool EncapsulatedGetdp::sol(const std::string options, const std::string modelName) {
//   onelab::server::citer it= onelab::server::instance()->findClient(getName());
//   onelab::client *c = it->second;
//   std::string commandLine = modelName + " " + options  ;
//   appendResolution(commandLine);
//   appendOption(commandLine,"-v",verbosity) ;
//   c->run(commandLine);
//   return true;
// }
// bool EncapsulatedGetdp::pos(const std::string options, const std::string modelName) {
//   onelab::server::citer it= onelab::server::instance()->findClient(getName());
//   onelab::client *c = it->second;
//   std::string commandLine = modelName + " " + options;
//   appendPostpro(commandLine);
//   appendOption(commandLine,"-v",verbosity) ;
//   c->run(commandLine);
//   return true;
// }

bool MetaModel::run(const std::string options, const std::string modelName) {
  onelab::server::citer it= onelab::server::instance()->findClient(getName());
  onelab::client *c = it->second;
  std::string commandLine = options  ;
  c->run(commandLine);
  return true;
}

bool MetaModel::analyze(const std::string options, const std::string modelName) {
  onelab::server::citer it= onelab::server::instance()->findClient(getName());
  onelab::client *c = it->second;
  std::string commandLine = "-a" + options  ;
  c->run(commandLine);
  return true;
}

/*
ONELAB additional tools
 */

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


#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
std::string getUserHomedir(){
  struct passwd *pw = getpwuid(getuid());
  std::string str(pw->pw_dir);
  str.append("/");
  return str;
}

#include <sys/param.h>
std::string getCurrentWorkdir(){
  char path[MAXPATHLEN];
  getcwd(path, MAXPATHLEN);
  std::string str = path;
  return str;
}

#include <sys/stat.h>		
#include <ctime>
bool checkIfPresent(std::string filename){
  struct stat buf;
  if (!stat(filename.c_str(), &buf))
    return true;
  else{
    Msg::Fatal("The file %s is not present",filename.c_str());
    return false;
  }
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
  printf("ONELAB System call(%s)\n", cmd.c_str());
  return system(cmd.c_str()); 
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
  if ( lin>=0 && lin<data.size())
    if (  col>=0 && col<data[lin-1].size())
      return data[lin-1][col-1];
  Msg::Fatal("The value has not been calculated: (%d,%d) out of range",lin,col);
}

array read_array(std::string filename, char sep){
  std::ifstream infile(filename.c_str());
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
    std::cerr << "Fooey!\n";
  }
  return array;
}


int getOptions(int argc, char *argv[],int &modelNumber, bool &analyzeOnly, std::string &sockName){
  int i = 1;
  while(i < argc) {
    if(argv[i][0] == '-') {
      if(!strcmp(argv[i] + 1, "m")) {
	i++;
	modelNumber = (atoi(argv[i]));
        i++;        
      }
      else if(!strcmp(argv[i] + 1, "onelab")) {
	i++;
	sockName = argv[i];
        i++;
      }
      else if(!strcmp(argv[i] + 1, "c")) {
	i++;
	std::cout << argv[0] << " has " << onelab::server::instance()->nbClients() << " clients\n" ;
	for(onelab::server::citer it = onelab::server::instance()->firstClient();
	    it != onelab::server::instance()->lastClient(); it++){
	  std::cout << it->second->getId() << ':' << it->second->getName() << std::endl;
	}
	exit(1);
      }
      else if(!strcmp(argv[i] + 1, "a")) {
	i++;
	analyzeOnly=true;
      }
      else {
	printf("Usage: %s [-m num -sol -c]\n", argv[0]);
	printf("Options are:\nm      model number\n");
	printf("c      list of clients\n");
	printf("a      analyze only\n");
	exit(1);
      }
    }
    else{
      std::string dummy=argv[i];
      i++;
    }
  }
}


  
