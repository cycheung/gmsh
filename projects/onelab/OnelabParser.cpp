#include "OnelabClients.h"
#include <algorithm>

std::string localSolverClient::longName(const std::string name){
  std::set<std::string>::iterator it;
  if((it = _parameters.find(name)) != _parameters.end())
    return *it;
  else
    return name;
}

std::string localSolverClient::resolveGetVal(std::string line) {
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
    get(strings,longName(arguments[0]));
    if (strings.size())
      buff.assign(strings[0].getValue());
    else
      Msg::Fatal("ONELAB unknown variable: %s",arguments[0].c_str());
    line.replace(pos0,pos-pos0,buff); 
    cursor=pos0+buff.length();
  }
  return line;
}

std::string localSolverClient::evaluateGetVal(std::string line) {
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  std::vector<std::string> arguments;
  std::string buff;
  int pos,cursor;

  if((pos=line.find(olkey::getValue,0)) == std::string::npos)
    return line;
  cursor = pos + olkey::getValue.length();
  pos=line.find(')',cursor);
  //std::cout << "ici:<" << line.substr(cursor,pos+1-cursor) << ">" << pos << std::endl;
  if(enclosed(line.substr(cursor,pos+1-cursor),arguments)<1)
    Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::getValue.c_str(),line.c_str());
  std::string paramName;
  paramName.assign(longName(arguments[0])); 
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
  return buff;
}


void localSolverClient::parse_oneline(std::string line, std::ifstream &infile) { 
  int pos0,pos,cursor;
  std::string name,action, path;
  std::vector<onelab::number> numbers;
  std::vector<std::string> arguments;
  std::vector<onelab::string> strings;
  char sep=';';
  std::string buff;
  std::set<std::string>::iterator it;

  if( (pos=line.find(olkey::client)) != std::string::npos) {// onelab.client
    parse_clientline(line,infile);
  }
  else if ( (pos=line.find(olkey::param)) != std::string::npos) {// onelab.param
    cursor = pos+olkey::param.length();
    while ( (pos=line.find(sep,cursor)) != std::string::npos){
      std::string name, action;
      extract(line.substr(cursor,pos-cursor),name,action,arguments);

      if(!action.compare("number")) { // paramName.number(val,path,help,...)
	if(arguments.empty())
	  Msg::Fatal("ONELAB: No value given for param <%s>",name.c_str());
	double val=atof(arguments[0].c_str());
	if(arguments.size()>1){
	  name.assign(arguments[1] + name);
	}
	_parameters.insert(name);
	get(numbers,name);
	if(numbers.empty()){ // creates parameter or skip if it exists
	  numbers.resize(1);
	  numbers[0].setName(name);
	  numbers[0].setValue(val);
	  if(arguments.size()>2)
	    numbers[0].setHelp(arguments[2]);
	  set(numbers[0]);
	}
      }
      else if(!action.compare("string")) { // paramName.string(val,path,help)
	if(arguments.empty())
	  Msg::Fatal("ONELAB: No value given for param <%s>",name.c_str());
	std::string value=arguments[0];
	if(arguments.size()>1)
	  name.assign(arguments[1] + name); // append path
	_parameters.insert(name);
	get(strings,name); 
	if(strings.empty()){ // creates parameter or skip if it exists
	  strings.resize(1);
	  strings[0].setName(name);
	  strings[0].setValue(value);
	  if(arguments.size()>2)
	    strings[0].setHelp(arguments[2]);
	  set(strings[0]);
	}
      }
      else if(!action.compare("AddChoices")){
	if(arguments.size()){
	  name.assign(longName(name));
	  get(numbers,name);
	  if(numbers.size()){ // parameter must exist
	    std::vector<double> choices=numbers[0].getChoices();
	    //numbers[0].setValue(atof(arguments[0].c_str()));
	    for(unsigned int i = 0; i < arguments.size(); i++){
	      double val=atof(arguments[i].c_str());
	      if(std::find(choices.begin(),choices.end(),val)==choices.end())
		choices.push_back(val);
	    }
	    numbers[0].setChoices(choices);
	    set(numbers[0]);
	  }
	  else{
	    get(strings,name);
	    if(strings.size()){
	      std::vector<std::string> choices=strings[0].getChoices();
	      //strings[0].setValue(arguments[0]);
	      for(unsigned int i = 0; i < arguments.size(); i++)
		if(std::find(choices.begin(),choices.end(),arguments[i])==choices.end())
		  choices.push_back(arguments[i]);
	      strings[0].setChoices(choices);
	      set(strings[0]);
	    }
	    else{
	      Msg::Fatal("ONELAB: the parameter <%s> does not exist",name.c_str());
	    }
	  }
	}
      }
      else if(!action.compare("MinMax")){
	if(arguments.size()){
	  name.assign(longName(name));
	  get(numbers,name);
	  if(numbers.size()){
	    numbers[0].setMin(atof(arguments[0].c_str()));
	    if(arguments.size()>1)
	      numbers[0].setMax(atof(arguments[1].c_str()));
	    if(arguments.size()>2)
	      numbers[0].setStep(atof(arguments[2].c_str()));
	    set(numbers[0]);
	  }
	}
      }
      else if(!action.compare("SetValue")){
	if(arguments.empty())
	  Msg::Fatal("ONELAB: missing argument SetValue <%s>",name.c_str());
	get(numbers,longName(name)); 
	if(numbers.size()){ 
	  numbers[0].setValue(atof(evaluateGetVal(arguments[0]).c_str()));
	  set(numbers[0]);
	}
	else{
	  get(strings,name); 
	  if(strings.size()){
	    strings[0].setValue(arguments[0]);
	    set(strings[0]);
	  }
	  else{
	    Msg::Fatal("ONELAB: the parameter <%s> does not exist",name.c_str());
	  }
	}
      }
      else
	Msg::Fatal("ONELAB: unknown action <%s>",action.c_str());
      cursor=pos+1;
    }
    // fin de la boucle while
  }
  else if ( (pos=line.find(olkey::iftrue)) != std::string::npos) {// onelab.iftrue
    cursor = pos+olkey::iftrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::iftrue.c_str(),line.c_str());
    bool condition = false;
    get(numbers,longName(arguments[0]));
    if (numbers.size())
      condition = (bool) numbers[0].getValue();
    if (!parse_ifstatement(infile,condition))
      Msg::Fatal("ONELAB misformed <%s> statement: <%s>",olkey::iftrue.c_str(),arguments[0].c_str());
  }
  else if ( (pos=line.find(olkey::ifequal)) != std::string::npos) {// onelab.ifequal
    cursor = pos+olkey::ifequal.length();
    pos=line.find_first_of(')',cursor)+1;
    if (enclosed(line.substr(cursor,pos-cursor),arguments) <2)
      Msg::Fatal("ONELAB misformed %s statement: <%s>",olkey::ifequal.c_str(),line.c_str());
    bool condition= false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= !strings[0].getValue().compare(arguments[1]);
    if (!parse_ifstatement(infile,condition))
      Msg::Fatal("ONELAB misformed <%s> statement: (%s,%s)",olkey::ifequal.c_str(),arguments[0].c_str(),arguments[1].c_str());
  }
  else if ( (pos=line.find(olkey::include)) != std::string::npos) {// onelab.include
    cursor = pos+olkey::include.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::include.c_str(),line.c_str());
    parse_onefile(arguments[0]);
  }
}

void localSolverClient::parse_onefile(std::string fileName) { 
  std::string line;
  std::string fileNameSave = fileName+".save";

  std::ifstream infile(fileName.c_str()); // read client description
  if (infile.is_open()){
    while ( infile.good() ) {
      getline (infile,line);
      parse_oneline(line,infile);
    }
    infile.close();
  }
  else
    Msg::Fatal("The file %s cannot be opened",fileName.c_str());
  infile.open(fileNameSave.c_str()); // read saved client pathes (if file present)
  if (infile.is_open()){
    while (infile.good() ) {
      getline (infile,line);
      parse_oneline(line,infile);
    }
    infile.close();
  }
} 


bool localSolverClient::parse_ifstatement(std::ifstream &infile, bool condition) { 
  int pos;
  std::string line;

  bool trueclause=true, statementend=false;
  while ( infile.good() && !statementend) {
    getline (infile,line);
    if ( (pos=line.find(olkey::olelse)) != std::string::npos) 
      trueclause=false;
    else if ( (pos=line.find(olkey::olendif)) != std::string::npos) 
      statementend=true;
    else if ( !(trueclause ^ condition) ) // xor bitwise operator
      parse_oneline(line,infile);
  }
  return statementend;
} 

bool localSolverClient::convert_ifstatement(std::ifstream &infile, std::ofstream &outfile, bool condition) { 
  int pos;
  std::string line;

  bool trueclause=true, statementend=false;
  while ( infile.good() && !statementend) {
    getline (infile,line);
    if ( (pos=line.find(olkey::olelse)) != std::string::npos) 
      trueclause=false;
    else if ( (pos=line.find(olkey::olendif)) != std::string::npos) 
      statementend=true;
    else if ( !(trueclause ^ condition) ) // xor bitwise operator
      convert_oneline(line,infile,outfile);
  }
  return statementend;
} 

void localSolverClient::convert_oneline(std::string line, std::ifstream &infile, std::ofstream &outfile) { 
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  std::vector<std::string> arguments;
  int pos0,pos,cursor;
  std::string buff;
  std::set<std::string>::iterator it;

  if ( (pos=line.find(olkey::label)) == std::string::npos) // not a onelab line
    outfile << line << std::endl; 
  else{ 
    if ( (pos=line.find(olkey::param)) != std::string::npos) {// onelab.param
      //skip the line
    }
    else if ( (pos=line.find(olkey::include)) != std::string::npos) {// onelab.include
      cursor = pos+olkey::include.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::include.c_str(),line.c_str());
      convert_onefile(arguments[0],outfile);
    }
    else if ( (pos=line.find(olkey::iftrue)) != std::string::npos) {// onelab.iftrue
      cursor = pos+olkey::iftrue.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::iftrue.c_str(),line.c_str());
      bool condition = false; 
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      if (!convert_ifstatement(infile,outfile,condition))
	Msg::Fatal("ONELAB misformed <%s> statement: %s",olkey::iftrue.c_str(),arguments[0].c_str());     
    }
    else if ( (pos=line.find(olkey::ifequal)) != std::string::npos) {// onelab.ifequal
      cursor = pos+olkey::ifequal.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<2)
	Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::ifequal.c_str(),line.c_str());;
      bool condition= false;
      get(strings,longName(arguments[0]));
      if (strings.size())
	condition =  !strings[0].getValue().compare(arguments[1]);
      if (!convert_ifstatement(infile,outfile,condition))
	Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::ifequal.c_str(),line.c_str());
    }
    else if ( (pos=line.find(olkey::getValue)) != std::string::npos) {// onelab.getValue
      // onelab.getValue, possibly several times in the line
      cursor=0;
      while ( (pos=line.find(olkey::getValue,cursor)) != std::string::npos){
	pos0=pos; // for further use
	cursor = pos+olkey::getValue.length();
	pos=line.find_first_of(')',cursor)+1;
	if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	  Msg::Fatal("ONELAB misformed <%s> statement: (%s)",olkey::getValue.c_str(),line.c_str());
	std::string paramName;
	paramName.assign(longName(arguments[0])); 
	//std::cout << "getValue:<" << arguments[0] << "> => <" << paramName << ">" << std::endl;
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
    else{
      outfile << line << std::endl; 
      Msg::Warning("ONELAB ambiguous sentence: %s",line.c_str());
    }
  }
}

void localSolverClient::convert_onefile(std::string ifilename, std::ofstream &outfile) { 
  std::string line;
  std::ifstream infile(ifilename.c_str());
  //fileName.assign(name.substr(0,name.find_last_of(".")));  // remove extension
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

void MetaModel::parse_clientline(std::string line, std::ifstream &infile) { 
  int pos,cursor;
  std::string name,action,path="";
  std::vector<std::string> arguments;
  std::vector<onelab::string> strings;
  char sep=';';

  if( (pos=line.find(olkey::client)) != std::string::npos) {// onelab.client
    cursor = pos + olkey::client.length();
    while ( (pos=line.find(sep,cursor)) != std::string::npos){
      extract(line.substr(cursor,pos-cursor),name,action,arguments);
      // 1: type, 2: path (optional)
      if(!action.compare("Register")){
	if(!findClientByName(name)){
	  Msg::Info("ONELAB: define client <%s>", name.c_str());
	  if(arguments.size()==1)
	    path="";
	  else if(arguments.size()==2)
	    path=arguments[1];
	  else
	    Msg::Error("ONELAB: wrong client definition <%s>", name.c_str());

	  if(path.empty()){ //check if one has a path on the server
	    get(strings,name + "/Path");
	    if(strings.size())
	      path=strings[0].getValue();
	  }
	  onelab::string o(name + "/Path",path);
	  o.setKind("file");
	  o.setVisible(path.empty());
	  o.setAttribute("Highlight","true");
	  set(o);
	  //client can be registered with path empty, but won't be run.
	  registerClient(name,arguments[0],path);
	}
	else
	  Msg::Error("ONELAB: redefinition of client <%s>", name.c_str());
      }
      else if(!action.compare("Path")){
	if(findClientByName(name)){
	  if(arguments.size()) {
	    if(arguments[0].size()){
	      onelab::string o(name + "/Path",arguments[0]);
	      o.setKind("file");
	      o.setVisible(false);
	      set(o);
	    }
	    else
	      Msg::Error("ONELAB: no path given for client <%s>", name.c_str());
	  }
	}
	else
	  Msg::Error("ONELAB: unknown client <%s>", name.c_str());
      }
      else if(!action.compare("Active")){
	localSolverClient *c;
	if(c=findClientByName(name)){
	  if(arguments.size()) {
	    if(arguments[0].size())
	      c->setActive(atof(arguments[0].c_str()));
	    else
	      Msg::Error("ONELAB: no path given for client <%s>", name.c_str());
	  }
	}
	else
	  Msg::Fatal("ONELAB: unknown client <%s>", name.c_str());
      }
      else if(!action.compare("Set")){
	if(arguments[0].size()){
	  strings.resize(1);
	  strings[0].setName(name);
	  strings[0].setValue(resolveGetVal(arguments[0]));
	  strings[0].setVisible(false);
	  if( (arguments[0].find(".geo") != std::string::npos) || 
              (arguments[0].find(".sif") != std::string::npos) ||
	      (arguments[0].find(".pro") != std::string::npos)) {
	    strings[0].setKind("file");
	    strings[0].setVisible(true);
	  }
	  std::vector<std::string> choices;
	  for(unsigned int i = 0; i < arguments.size(); i++)
	    if(std::find(choices.begin(),choices.end(),arguments[i])==choices.end())
	      choices.push_back(resolveGetVal(arguments[i]));
	  strings[0].setChoices(choices);
	  set(strings[0]);
	}
      }
      else if(!action.compare("List")){//no check whether choices[i] already inserted
	if(arguments[0].size()){
	  strings.resize(1);
	  strings[0].setName(name);
	  strings[0].setValue(resolveGetVal(arguments[0]));
	  strings[0].setVisible(false);
	  std::vector<std::string> choices;
	  for(unsigned int i = 0; i < arguments.size(); i++)
	    choices.push_back(resolveGetVal(arguments[i]));
	  strings[0].setChoices(choices);
	  set(strings[0]);
	}
      }
      else
	Msg::Fatal("ONELAB:unknown keyword <%s>",action.c_str());
      cursor=pos+1;
    }
  }
}
