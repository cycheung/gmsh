#include "OnelabClients.h"
#include <algorithm>

// reserved keywords for the onelab parser

namespace olkey{ 
  static std::string deflabel("onelab.tags");
  static std::string label("OL."), comment("%"), separator(";");
  static std::string client(label+"client");
  static std::string param(label+"parameter");
  static std::string number(label+"number"), string(label+"string");
  static std::string include(label+"include");
  static std::string iftrue(label+"iftrue"), olelse(label+"else"), olendif(label+"endif");
  static std::string ifequal(label+"ifequal"), ifntrue(label+"ifntrue");
  static std::string getValue(label+"getValue");
  static std::string arguments("Args"), inFiles("In"), outFiles("Out");
}

std::string localSolverClient::toChar(){
  std::ostringstream sstream;
  if(getCommandLine().size()){
    sstream << olkey::client << " " 
	    << getName() << "." << "CommandLine(" 
	    << getCommandLine() << ");\n";
  }
  return sstream.str();
}

std::string localSolverClient::longName(const std::string name){
  std::set<std::string>::iterator it;
  if((it = _parameters.find(name)) != _parameters.end())
    return *it;
  else
    return name;
}

std::string localSolverClient::resolveGetVal(std::string line) {
  std::vector<onelab::number> numbers;
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
      Msg::Fatal("Misformed <%s> statement: (%s)",olkey::getValue.c_str(),line.c_str());
    std::string paramName=longName(arguments[0]);
    get(numbers,paramName);
    if (numbers.size()){
      std::stringstream Num;
      Num << numbers[0].getValue();
      buff.assign(Num.str());
    }
    else{
      get(strings,longName(paramName));
      if (strings.size())
	buff.assign(strings[0].getValue());
      else
	Msg::Fatal("resolveGetVal: unknown variable: %s",arguments[0].c_str());
    }
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
    Msg::Fatal("Misformed <%s> statement: (%s)",olkey::getValue.c_str(),line.c_str());
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
      Msg::Fatal("evaluetaGetVal: unknown variable: %s",paramName.c_str());
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

  if((pos=line.find_first_not_of(" "))==std::string::npos)
    return; // empty line
  else if(!line.compare(pos,olkey::comment.size(),olkey::comment))
    return; // commented out, skip the line
  else if ( (pos=line.find(olkey::deflabel)) != std::string::npos){
    // onelab.tags(label,comment,separator)
    cursor = pos+olkey::deflabel.length();
    pos=line.find_first_of(')',cursor)+1;
    int NumArg=enclosed(line.substr(cursor,pos-cursor),arguments);
    if (NumArg==0)
      Msg::Fatal("No tag given, using default tags <%s,%s,%s>",
		 olkey::label.c_str(),olkey::comment.c_str()),olkey::separator.c_str();
    else if(NumArg>=1){
      olkey::label.assign(arguments[0]);
      olkey::client.assign(olkey::label+"client");
      olkey::param.assign(olkey::label+"parameter");
      olkey::number.assign(olkey::label+"number"), olkey::string.assign(olkey::label+"string");
      olkey::include.assign(olkey::label+"include"); 
      olkey::iftrue.assign(olkey::label+"iftrue"), olkey::olelse.assign(olkey::label+"else");
      olkey::olendif.assign(olkey::label+"endif"); 
      olkey::ifequal.assign(olkey::label+"ifequal");
      olkey::getValue.assign(olkey::label+"getValue");
      std::cout << "Tags modified:" << olkey::client << std::endl;
    }
    else if(NumArg>=2)
      olkey::comment.assign(arguments[1]);
    else if(NumArg>=3)
      olkey::separator.assign(arguments[2]);
  }
  else if ( (pos=line.find(olkey::label)) == std::string::npos) // not a onelab line
    return;
  else if( (pos=line.find(olkey::client)) != std::string::npos) {// onelab.client
    parse_clientline(line,infile);
  }
  else if ( (pos=line.find(olkey::param)) != std::string::npos) {// onelab.param
    cursor = pos+olkey::param.length();
    while ( (pos=line.find(sep,cursor)) != std::string::npos){
      std::string name, action;
      extract(line.substr(cursor,pos-cursor),name,action,arguments);

      if(!action.compare("number")) { 
	// syntax: paramName.number(val,path,help,range(optional))
	if(arguments.empty())
	  Msg::Fatal("No value given for param <%s>",name.c_str());
	double val=atof(arguments[0].c_str());
	if(arguments.size()>=2){
	  name.assign(arguments[1] + name);
	}
	_parameters.insert(name);
	get(numbers,name);
	if(numbers.empty()) { //if already exists, skip 
	  numbers.resize(1);
	  numbers[0].setName(name);
	  numbers[0].setValue(val);
	  if(arguments.size()>=3)
	    numbers[0].setHelp(arguments[2]);
	  if(arguments.size()>=4){
	    std::vector<double> bounds;
	    if (extractRange(arguments[3],bounds)){
	      numbers[0].setMin(bounds[0]);
	      numbers[0].setMax(bounds[1]);
	      numbers[0].setStep(bounds[2]);
	    }
	  }
	  set(numbers[0]);
	}
      }
      else if(!action.compare("MinMax")){ // add a range to an existing prameter
	if(arguments.empty())
	  Msg::Fatal("No argument given for MinMax <%s>",name.c_str());
	name.assign(longName(name));
	get(numbers,name);
	bool noRange = true;
	if(numbers.size()){ //parameter must exist
	  // check whether a range already exists
	  if(numbers[0].getMin() != -onelab::parameter::maxNumber() ||
	     numbers[0].getMax() != onelab::parameter::maxNumber() ||
	     numbers[0].getStep() != 0.) noRange = false;
	  if(noRange){ //if a range is already defined, skip
	    if(arguments.size()==1){
	      std::vector<double> bounds;
	      if (extractRange(arguments[1],bounds)){
		numbers[0].setMin(bounds[0]);
		numbers[0].setMax(bounds[1]);
		numbers[0].setStep(bounds[2]);
	      }
	    }
	    else if(arguments.size()==3){
	      numbers[0].setMin(atof(arguments[0].c_str()));
	      numbers[0].setMax(atof(arguments[1].c_str()));
	      numbers[0].setStep(atof(arguments[2].c_str()));
	    }
	    else
	      Msg::Fatal("Wrong argument number for MinMax <%s>",name.c_str());
	    set(numbers[0]);
	  }
	}
      }
      else if(!action.compare("string")) { // paramName.string(val,path,help)
	if(arguments.empty())
	  Msg::Fatal("No value given for param <%s>",name.c_str());
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
	      Msg::Fatal("The parameter <%s> does not exist",name.c_str());
	    }
	  }
	}
      }
      else if(!action.compare("SetValue")){
	if(arguments.empty())
	  Msg::Fatal("Missing argument SetValue <%s>",name.c_str());
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
	    Msg::Fatal("The parameter <%s> does not exist",name.c_str());
	  }
	}
      }
      else
	Msg::Fatal("Unknown action <%s>",action.c_str());
      cursor=pos+1;
    }
    // fin de la boucle while
  }
  else if ( (pos=line.find(olkey::iftrue)) != std::string::npos) {// onelab.iftrue
    cursor = pos+olkey::iftrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",olkey::iftrue.c_str(),line.c_str());
    bool condition = false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= !strings[0].getValue().compare(arguments[1]);
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      if (!parse_ifstatement(infile,condition))
	Msg::Fatal("Misformed <%s> statement: <%s>",olkey::iftrue.c_str(),arguments[0].c_str());
    }
  }
  else if ( (pos=line.find(olkey::ifntrue)) != std::string::npos) {// onelab.ifntrue
    cursor = pos+olkey::ifntrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",olkey::ifntrue.c_str(),line.c_str());
    bool condition = false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= !strings[0].getValue().compare(arguments[1]);
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      if (!parse_ifstatement(infile,!condition))
	Msg::Fatal("Misformed <%s> statement: <%s>",olkey::iftrue.c_str(),arguments[0].c_str());
    }
  }
  else if ( (pos=line.find(olkey::ifequal)) != std::string::npos) {// onelab.ifequal
    cursor = pos+olkey::ifequal.length();
    pos=line.find_first_of(')',cursor)+1;
    if (enclosed(line.substr(cursor,pos-cursor),arguments) <2)
      Msg::Fatal("Misformed %s statement: <%s>",olkey::ifequal.c_str(),line.c_str());
    bool condition= false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= !strings[0].getValue().compare(arguments[1]);
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition= (numbers[0].getValue() == atof(arguments[1].c_str()));
      else
	Msg::Fatal("Unknown argument <%s> in <%s> statement", arguments[0].c_str()),olkey::ifequal.c_str();
    }
    if (!parse_ifstatement(infile,condition))
      Msg::Fatal("Misformed <%s> statement: (%s,%s)",
	         olkey::ifequal.c_str(),arguments[0].c_str(),arguments[1].c_str());
  }
  else if ( (pos=line.find(olkey::include)) != std::string::npos) { // onelab.include
    cursor = pos+olkey::include.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",olkey::include.c_str(),line.c_str());
    parse_onefile(arguments[0]);
  }
}

void localSolverClient::parse_onefile(std::string fileName) { 
  std::string line;
  // std::string fileNameSave = fileName+".save";

  std::ifstream infile(fileName.c_str()); // read client description
  if (infile.is_open()){
    Msg::Info("Read file <%s>",fileName.c_str());
    while ( infile.good() ) {
      getline (infile,line);
      parse_oneline(line,infile);
    }
    infile.close();
  }
  else
    Msg::Info("The file %s does not exist",fileName.c_str());
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

  if((pos=line.find_first_not_of(" "))==std::string::npos)
    outfile << line << std::endl;  // empty line
  else if(!line.compare(pos,olkey::comment.size(),olkey::comment))
    return;  // commented out, skip the line
  else if ( (pos=line.find(olkey::deflabel)) != std::string::npos)
    return;
  else if ( (pos=line.find(olkey::label)) == std::string::npos) 
    outfile << line << std::endl; // not a onelab line
  else{ 
    if(!line.compare(line.find_first_not_of(" "),olkey::comment.size(),olkey::comment)){
      // commented out, skip the line
    }
    else if ( (pos=line.find(olkey::param)) != std::string::npos) {// onelab.param
      //skip the line
    }
    else if ( (pos=line.find(olkey::include)) != std::string::npos) {// onelab.include
      cursor = pos+olkey::include.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	Msg::Fatal("Misformed <%s> statement: (%s)",olkey::include.c_str(),line.c_str());
      convert_onefile(arguments[0],outfile);
    }
    else if ( (pos=line.find(olkey::iftrue)) != std::string::npos) {// onelab.iftrue
      cursor = pos+olkey::iftrue.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	Msg::Fatal("Misformed <%s> statement: (%s)",olkey::iftrue.c_str(),line.c_str());
      bool condition = false; 
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      if (!convert_ifstatement(infile,outfile,condition))
	Msg::Fatal("Misformed <%s> statement: %s",olkey::iftrue.c_str(),arguments[0].c_str());     
    }
    else if ( (pos=line.find(olkey::ifequal)) != std::string::npos) {// onelab.ifequal
      cursor = pos+olkey::ifequal.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<2)
	Msg::Fatal("Misformed <%s> statement: (%s)",olkey::ifequal.c_str(),line.c_str());;
      bool condition= false;
      get(strings,longName(arguments[0]));
      if (strings.size())
	condition =  !strings[0].getValue().compare(arguments[1]);
      if (!convert_ifstatement(infile,outfile,condition))
	Msg::Fatal("Misformed <%s> statement: (%s)",olkey::ifequal.c_str(),line.c_str());
    }
    else if ( (pos=line.find(olkey::getValue)) != std::string::npos) {// onelab.getValue
      // onelab.getValue, possibly several times in the line
      cursor=0;
      while ( (pos=line.find(olkey::getValue,cursor)) != std::string::npos){
	pos0=pos; // for further use
	cursor = pos+olkey::getValue.length();
	pos=line.find_first_of(')',cursor)+1;
	if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	  Msg::Fatal("Misformed <%s> statement: (%s)",olkey::getValue.c_str(),line.c_str());
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
	    Msg::Fatal("Unknown variable: %s",paramName.c_str());
	}
	line.replace(pos0,pos-pos0,buff); 
	cursor=pos0+buff.length();
      }
      outfile << line << std::endl; 
    }
    else{
      outfile << line << std::endl; 
      Msg::Warning("Ambiguous sentence: %s",line.c_str());
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
  std::string name,action,cmdl;
  std::vector<std::string> arguments;
  std::vector<onelab::string> strings;
  char sep=';';

  if( (pos=line.find(olkey::client)) != std::string::npos) {// onelab.client
    // syntax name.Register(interfaced)
    //        name.Register(interfaced,cmdl)
    //        name.Register(interfaced,remoteHost,remoteDir,cmdl)
    //        name.Register(encapsulated)
    //        name.Register(encapsulated,remoteHost,remoteDir,cmdl)
    cursor = pos + olkey::client.length();
    while ( (pos=line.find(sep,cursor)) != std::string::npos){
      extract(line.substr(cursor,pos-cursor),name,action,arguments);
      if(!action.compare("Register")){
	if(!findClientByName(name)){
	  Msg::Info("Define client <%s>", name.c_str());

	  if(arguments.size()==1){ //local clients
	    //check if one has a saved command line on the server (prealably read from file .ol.save)
	    get(strings,name + "/CommandLine");
	    if(strings.size())
	      cmdl.assign(strings[0].getValue());
	    else 
	      cmdl.assign("");
	    registerClient(name,resolveGetVal(arguments[0]),cmdl); // possibly with an empty cmdl
	  }
	  else if(arguments.size()==2){ //local clients
	    registerClient(name,resolveGetVal(arguments[0]),resolveGetVal(arguments[1])); 
	  }
	  else if(arguments.size()==4){ //remote clients, disregard saved commandLine
	    cmdl.assign(resolveGetVal(arguments[3]));
	    registerClient(name,resolveGetVal(arguments[0]),cmdl,
			   resolveGetVal(arguments[1]),resolveGetVal(arguments[2]));
	  }
	  else
	    Msg::Fatal("Wrong number or arguments in client definition <%s>", name.c_str());

	  //set actual commandLine on server (for consistency)
	  onelab::string o(name + "/CommandLine","");
	  o.setValue(cmdl);
	  o.setKind("file");
	  o.setVisible(cmdl.empty());
	  o.setAttribute("Highlight","true");
	  set(o);
	}
	else
	  Msg::Error("Redefinition of client <%s>", name.c_str());
      }
      else if(!action.compare("CommandLine")){
	if(findClientByName(name)){
	  if(arguments.size()) {
	    if(arguments[0].size()){
	      onelab::string o(name + "/CommandLine",arguments[0]);
	      o.setKind("file");
	      o.setVisible(false);
	      set(o);
	    }
	    else
	      Msg::Error("No pathname given for client <%s>", name.c_str());
	  }
	}
	else
	  Msg::Error("Unknown client <%s>", name.c_str());
      }
      else if(!action.compare("Active")){
	localSolverClient *c;
	if(c=findClientByName(name)){
	  if(arguments.size()) {
	    if(arguments[0].size())
	      c->setActive(atof(arguments[0].c_str()));
	    else
	      Msg::Error("No argument for <%s.Active> statement", name.c_str());
	  }
	}
	else
	  Msg::Fatal("Unknown client <%s>", name.c_str());
      }
      else if(!action.compare("Args")){
	if(arguments[0].size()){
	  strings.resize(1);
	  strings[0].setName(name+"/Arguments");
	  strings[0].setValue(resolveGetVal(arguments[0]));
	  strings[0].setVisible(false);
	  set(strings[0]);
	}
      }
      else if(!action.compare("In")){
	if(arguments[0].size()){
	  strings.resize(1);
	  strings[0].setName(name+"/InputFiles");
	  strings[0].setValue(resolveGetVal(arguments[0]));
	  strings[0].setKind("file");
	  strings[0].setVisible(false);
	  std::vector<std::string> choices;
	  for(unsigned int i = 0; i < arguments.size(); i++)
	    //if(std::find(choices.begin(),choices.end(),arguments[i])==choices.end())
	    choices.push_back(resolveGetVal(arguments[i]));
	  strings[0].setChoices(choices);
	  set(strings[0]);
	}
      }
      else if(!action.compare("Out")){
	if(arguments[0].size()){
	  strings.resize(1);
	  strings[0].setName(name+"/OutputFiles");
	  strings[0].setValue(resolveGetVal(arguments[0]));
	  strings[0].setKind("file");
	  strings[0].setVisible(false);
	  std::vector<std::string> choices;
	  for(unsigned int i = 0; i < arguments.size(); i++)
	    //if(std::find(choices.begin(),choices.end(),arguments[i])==choices.end())
	    choices.push_back(resolveGetVal(arguments[i]));
	  strings[0].setChoices(choices);
	  set(strings[0]);
	}
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
	Msg::Fatal("Unknown keyword <%s>",action.c_str());
      cursor=pos+1;
    }
  }
}
