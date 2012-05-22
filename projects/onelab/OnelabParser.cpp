#include "OnelabClients.h"
#include <algorithm>

// reserved keywords for the onelab parser

namespace olkey{ 
  static std::string deflabel("onelab.tags");
  static std::string label("OL."), comment("%"), separator(";");
  static std::string line(label+"line");
  static std::string begin(label+"begin");
  static std::string end(label+"end");
  static std::string setValue(label+"setValue");
  static std::string include(label+"include");
  static std::string ifcond(label+"if");
  static std::string ifequal(label+"ifequal");
  static std::string iftrue(label+"iftrue"), ifntrue(label+"ifntrue");
  static std::string olelse(label+"else"), olendif(label+"endif");
  static std::string getValue(label+"getValue");
  static std::string getRegion(label+"getRegion");
  static std::string number("number"), string("string");
  static std::string arguments("Args"), inFiles("In"), outFiles("Out");
  static std::string upload("Up"), merge("Merge");
  static std::string checkCmd("Check"), computeCmd("Compute");
}

int enclosed(const std::string &in, std::vector<std::string> &arguments){
  // syntax: (arguments[Ø], arguments[1], ... , arguments[n])
  int pos, cursor;
  arguments.resize(0);
  cursor=0;
  if ( (pos=in.find("(",cursor)) == std::string::npos )
     Msg::Fatal("Syntax error: <%s>",in.c_str());

  unsigned int count=1;
  pos++; // skips '('
  cursor = pos; 
  do{
    if(in[pos]=='(') count++;
    if(in[pos]==')') count--;
    if(in[pos]==',') {
      arguments.push_back(removeBlanks(in.substr(cursor,pos-cursor)));
      if(count!=1)
	Msg::Fatal("Syntax error: <%s>",in.c_str());
      cursor=pos+1; // skips ','
    }
    pos++;
  } while( count && (pos!=std::string::npos) );
  // count is 0 when the closing brace is found. 

  if(count)
     Msg::Fatal("Syntax error: <%s>",in.c_str());
  else
    arguments.push_back(removeBlanks(in.substr(cursor,pos-1-cursor)));
  return arguments.size();
}

int extractLogic(const std::string &in, std::vector<std::string> &arguments){
  // syntax: ( argument[0], argument[1]\in{<,>,<=,>=,==,!=}, arguments[2])
  int pos, cursor;
  arguments.resize(0);
  cursor=0;
  if ( (pos=in.find("(",cursor)) == std::string::npos )
     Msg::Fatal("Syntax error: <%s>",in.c_str());

  unsigned int count=1;
  pos++; // skips '('
  cursor=pos; 
  do{
    if(in[pos]=='(') count++;
    if(in[pos]==')') count--;
    if( (in[pos]=='<') || (in[pos]=='=') || (in[pos]=='>') ){
      arguments.push_back(removeBlanks(in.substr(cursor,pos-cursor)));
      if(count!=1)
	Msg::Fatal("Syntax error: <%s>",in.c_str());
      cursor=pos;
      if(in[pos+1]=='='){
	arguments.push_back(in.substr(cursor,2));
	pos++;
      }
      else{
      	arguments.push_back(in.substr(cursor,1));
      }
      cursor=pos+1;
    }
    pos++;
  } while( count && (pos!=std::string::npos) );
  // count is 0 when the closing brace is found. 

  if(count)
     Msg::Fatal("Syntax error: mismatched parenthesis in <%s>",in.c_str());
  else
    arguments.push_back(removeBlanks(in.substr(cursor,pos-1-cursor)));

  if((arguments.size()!=1) && (arguments.size()!=3))
    Msg::Fatal("Syntax error: <%s>",in.c_str());
  return arguments.size();
}

int extract(const std::string &in, std::string &paramName, std::string &action, std::vector<std::string> &arguments){
  // syntax: paramName.action( arg1, arg2, ... )
  int pos, cursor,NumArg=0;
  cursor=0;
  if ( (pos=in.find(".",cursor)) == std::string::npos )
     Msg::Fatal("Syntax error: <%s>",in.c_str());
  else
    paramName.assign(sanitize(in.substr(cursor,pos-cursor)));
  cursor = pos+1; // skips '.'
  if ( (pos=in.find("(",cursor)) == std::string::npos )
     Msg::Fatal("Syntax error: <%s>",in.c_str());
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
     Msg::Fatal("Syntax error: %s",in.c_str());
  else
    NumArg = enclosed(in.substr(cursor,pos-cursor),arguments);
  // std::cout << "paramName=<" << paramName << ">" << std::endl;
  // std::cout << "arguments=<" << in.substr(cursor,pos+1-cursor) << ">" << std::endl;
  return NumArg;
}

bool extractRange(const std::string &in, std::vector<double> &arguments){
  // syntax: a:b:c or a:b#n
  int pos, cursor;
  arguments.resize(0);
  cursor=0;
  if ( (pos=in.find(":",cursor)) == std::string::npos )
     Msg::Fatal("Syntax error in range <%s>",in.c_str());
  else{
    arguments.push_back(atof(in.substr(cursor,pos-cursor).c_str()));
  }
  cursor = pos+1; // skips ':'
  if ( (pos=in.find(":",cursor)) != std::string::npos ){
    arguments.push_back(atof(in.substr(cursor,pos-cursor).c_str()));
    arguments.push_back(atof(in.substr(pos+1).c_str()));
  }
  else if ( (pos=in.find("#",cursor)) != std::string::npos ){
    arguments.push_back(atof(in.substr(cursor,pos-cursor).c_str()));
    double NumStep = atof(in.substr(pos+1).c_str());
    arguments.push_back((arguments[1]-arguments[0])/((NumStep==0)?1:NumStep));
  }
  else
     Msg::Fatal("Syntax error in range <%s>",in.c_str());
  return (arguments.size()==3);
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
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::getValue.c_str(),line.c_str());
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
	Msg::Fatal("resolveGetVal: unknown variable: %s",paramName.c_str());
    }
    line.replace(pos0,pos-pos0,buff); 
    cursor=pos0+buff.length();
  }
  return line; // appeler ici éventuellement mathex, strcmp
}

bool localSolverClient::resolveLogicExpr(std::vector<std::string> arguments) {
  std::vector<onelab::number> numbers;
  double val1, val2;
  bool condition;

  val1= atof( resolveGetVal(arguments[0]).c_str() );
  if(arguments.size()==1)
    return (bool)val1;

  if(arguments.size()==3){
    val2=atof( resolveGetVal(arguments[2]).c_str() );
    if(!arguments[1].compare("<"))
      condition = (val1<val2);
    else if (!arguments[1].compare("<="))
      condition = (val1<=val2);
    else if (!arguments[1].compare(">"))
      condition = (val1>val2);
    else if (!arguments[1].compare(">="))
      condition = (val1>=val2);
    else if (!arguments[1].compare("=="))
      condition = (val1==val2);   
    else if (!arguments[1].compare("!="))
      condition = (val1!=val2);
  }

  return condition;
}

void localSolverClient::parse_sentence(std::string line) { 
  int pos,cursor;
  std::string name,action,path;
  std::vector<std::string> arguments;
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;

  cursor = 0;
  while ( (pos=line.find(olkey::separator,cursor)) != std::string::npos){
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
	  numbers[0].setLabel(arguments[2]);
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
    else if(!action.compare("MinMax")){ 
      // add a range to an existing number
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
    else if(!action.compare("string")) { 
      // paramName.string(val,path,help)
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
	  strings[0].setLabel(arguments[2]);
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
	      if(std::find(choices.begin(),choices.end(),
			   arguments[i])==choices.end())
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
    else if(!action.compare("AddLabels")){
      if(arguments.size()){
	name.assign(longName(name));
	get(numbers,name);
	if(numbers.size()){ // parameter must exist
	  std::vector<double> choices=numbers[0].getChoices();
	  if(choices.size() != arguments.size())
	    Msg::Fatal("Nb of labels does not match nb of choices <%s>",
		       name.c_str());
	  std::vector<std::string> labels;
	  for(unsigned int i = 0; i < arguments.size(); i++){
	      labels.push_back(arguments[i]);
	  }
	  std::cout << "calling setChoicesLabel" << std::endl;
	  numbers[0].setChoiceLabels(labels);
	  set(numbers[0]);
	}
	else
	  Msg::Fatal("The number <%s> does not exist",name.c_str());
      }
    }
    else if(!action.compare("SetValue")){
      if(arguments.empty())
	Msg::Fatal("Missing argument SetValue <%s>",name.c_str());
      get(numbers,longName(name)); 
      if(numbers.size()){ 
	numbers[0].setValue(atof(resolveGetVal(arguments[0]).c_str()));
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
    else if(!action.compare("SetVisible")){
      if(arguments.empty())
	Msg::Fatal("Missing argument SetVisible <%s>",name.c_str());
      get(numbers,longName(name)); 
      if(numbers.size()){ 
	numbers[0].setVisible(atof(resolveGetVal(arguments[0]).c_str()));
	set(numbers[0]);
      }
      else{
	get(strings,name); 
	if(strings.size()){
	  strings[0].setVisible(atof(resolveGetVal(arguments[0]).c_str()));
	  set(strings[0]);
	}
	else{
	  Msg::Fatal("The parameter <%s> does not exist",name.c_str());
	}
      }
    }
    else if(!action.compare("SetReadOnly")){
      if(arguments.empty())
	Msg::Fatal("Missing argument SetReadOnly <%s>",name.c_str());
      get(numbers,longName(name)); 
      if(numbers.size()){ 
	numbers[0].setReadOnly(atof(resolveGetVal(arguments[0]).c_str()));
	set(numbers[0]);
      }
      else{
	get(strings,name); 
	if(strings.size()){
	  strings[0].setReadOnly(atof(resolveGetVal(arguments[0]).c_str()));
	  set(strings[0]);
	}
	else{
	  Msg::Fatal("The parameter <%s> does not exist",name.c_str());
	}
      }
    }
    else
      client_sentence(name,action,arguments);
    cursor=pos+1;
  }
}

void localSolverClient::parse_oneline(std::string line, std::ifstream &infile) { 
  int pos,cursor;
  std::vector<std::string> arguments;
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;

  if((pos=line.find_first_not_of(" \t"))==std::string::npos){
    // empty line, skip
  }
  else if(!line.compare(pos,olkey::comment.size(),olkey::comment)){
    // commented out line, skip
  }
  else if ( (pos=line.find(olkey::deflabel)) != std::string::npos){
    // onelab.tags(label,comment,separator)
    cursor = pos+olkey::deflabel.length();
    pos=line.find_first_of(')',cursor)+1;
    int NumArg=enclosed(line.substr(cursor,pos-cursor),arguments);
    if (NumArg == 0)
      Msg::Error("No onelab tags given");
    if(NumArg >= 1){
      olkey::label.assign(arguments[0]);
      olkey::line.assign(olkey::label+"line");
      olkey::begin.assign(olkey::label+"begin");
      olkey::end.assign(olkey::label+"end");
      olkey::include.assign(olkey::label+"include");
      olkey::ifcond.assign(olkey::label+"if");
      olkey::ifequal.assign(olkey::label+"ifequal");
      olkey::iftrue.assign(olkey::label+"iftrue");
      olkey::ifntrue.assign(olkey::label+"ifntrue"); 
      olkey::olelse.assign(olkey::label+"else");
      olkey::olendif.assign(olkey::label+"endif");
      olkey::getValue.assign(olkey::label+"getValue");
      olkey::getRegion.assign(olkey::label+"getRegion");
    }
    if(NumArg >= 2)
      olkey::comment.assign(arguments[1]);
    if(NumArg >= 3)
      olkey::separator.assign(arguments[2]);
    Msg::Info("Using now onelab tags <%s,%s,%s>",
	      olkey::label.c_str(), olkey::comment.c_str(),
	      olkey::separator.c_str());
  }
  else if( (pos=line.find(olkey::begin)) != std::string::npos) {
    // onelab.begin
    if (!parse_block(infile))
      Msg::Fatal("Misformed <%s> block <%s>",
		 olkey::begin.c_str(),olkey::end.c_str());
  }
  else if ( (pos=line.find(olkey::iftrue)) != std::string::npos) {
    // onelab.iftrue
    cursor = pos+olkey::iftrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::iftrue.c_str(),line.c_str());
    bool condition = false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= true;
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      else
	Msg::Fatal("Unknown parameter <%s> in <%s> statement",
		   arguments[0].c_str(),olkey::iftrue.c_str());
      if (!parse_ifstatement(infile,condition))
	Msg::Fatal("Misformed <%s> statement: <%s>",
		   olkey::iftrue.c_str(),arguments[0].c_str());
    }
  }
  else if ( (pos=line.find(olkey::ifntrue)) != std::string::npos) {
    // onelab.ifntrue
    cursor = pos+olkey::ifntrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::ifntrue.c_str(),line.c_str());
    bool condition = false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= true;
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      else
	Msg::Fatal("Unknown parameter <%s> in <%s> statement",
		   arguments[0].c_str(),olkey::ifntrue.c_str());
      if (!parse_ifstatement(infile,!condition))
	Msg::Fatal("Misformed <%s> statement: <%s>",
		   olkey::ifntrue.c_str(),arguments[0].c_str());
    }
  }
  else if ( (pos=line.find(olkey::ifcond)) != std::string::npos) {
    // onelab.ifcond
    cursor = pos+olkey::ifcond.length();
    int NumArgs=extractLogic(line.substr(cursor),arguments);
    bool condition= resolveLogicExpr(arguments);
    if (!parse_ifstatement(infile,condition))
      Msg::Fatal("Misformed %s statement: <%s>", line.c_str());
  }
  else if ( (pos=line.find(olkey::ifequal)) != std::string::npos) {
    // onelab.ifequal
    cursor = pos+olkey::ifequal.length();
    pos=line.find_first_of(')',cursor)+1;
    if (enclosed(line.substr(cursor,pos-cursor),arguments) <2)
      Msg::Fatal("Misformed %s statement: <%s>",
		 olkey::ifequal.c_str(),line.c_str());
    bool condition= false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= !strings[0].getValue().compare(arguments[1]);
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
        condition= (numbers[0].getValue() == atof(arguments[1].c_str()));
      else
        Msg::Fatal("Unknown argument <%s> in <%s> statement",
		   arguments[0].c_str()),olkey::ifequal.c_str();
    }
    if (!parse_ifstatement(infile,condition))
      Msg::Fatal("Misformed <%s> statement: (%s,%s)",
                 olkey::ifequal.c_str(), arguments[0].c_str(),
		 arguments[1].c_str());
  }
  else if ( (pos=line.find(olkey::include)) != std::string::npos) { 
    // onelab.include
    cursor = pos+olkey::include.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::include.c_str(),line.c_str());
    parse_onefile(arguments[0]);
  }
  else if( isOnelabBlock() ||
	 ( !isOnelabBlock() &&
	   ((pos=line.find(olkey::line)) != std::string::npos)) ){
    // either any other line within a onelabBlock or a line 
    // introduced by a "onelab.line" tag not within a onelabBlock
    std::string cmds="",cmd;
    int posa, posb, NbLines=1;
    do{
      if( (pos=line.find(olkey::line)) != std::string::npos)
	posa=pos + olkey::line.size();
      else
	posa=0; // skip tag 'olkey::line' if any

      posb=line.find(olkey::comment); // skip trailing comments if any
      if(posb==std::string::npos)
	cmd.assign(line.substr(posa));
      else
	cmd.assign(line.substr(posa,posb-posa));
      cmds.append(cmd);

      //std::cout << "cmds=<" << cmds << ">" << std::endl;

      // check whether "cmd" ends now with "olkey::separator"
      posb=cmd.find_last_not_of(" \t")-olkey::separator.length()+1;
      if(posb<0) posb=0;
      if(cmd.compare(posb,olkey::separator.length(),olkey::separator)){
	// append the next line
	getline (infile,line);
	if((pos=line.find_first_not_of(" \t"))==std::string::npos){
	  Msg::Fatal("Empty line not allowed within a command <%s>",
		     cmds.c_str());
	}
	else if(!line.compare(pos,olkey::comment.size(),olkey::comment)){
	  Msg::Fatal("Comment lines not allowed within a command <%s>",
		     cmds.c_str());
	}
	NbLines++; // command should not span over more than 10 lines
      }
      else
	break;
    } while (infile.good() && NbLines<=10);
    if(NbLines>=10)
      Msg::Fatal("Command <%s> should not span over more than 10 lines",
		 cmds.c_str());
    parse_sentence(cmds);
  }
  else if ( (pos=line.find(olkey::getValue)) != std::string::npos) {
    // onelab.getValue: nothing to do
  }
  else if ( (pos=line.find(olkey::getRegion)) != std::string::npos) {
    // onelab.getRegion: nothing to do
  }
  else if( (pos=line.find(olkey::label)) != std::string::npos) {
      Msg::Fatal("Unknown ONELAB keyword in <%s>",line.c_str());
  }
  else{
    // not a onelab line, skip
  }
}

bool localSolverClient::parse_block(std::ifstream  &infile) { 
  int level, pos;
  std::string line;
  openOnelabBlock();
  while (infile.good()){
    getline (infile,line);
    if ((pos=line.find(olkey::end)) != std::string::npos){
      closeOnelabBlock();
      return true;
    }
    parse_oneline(line,infile);
  }
  return false;
} 

bool localSolverClient::parse_ifstatement(std::ifstream &infile, 
					  bool condition) { 
  int level, pos;
  std::string line;

  bool trueclause=true; 
  level=1;
  while ( infile.good() && level) {
    getline (infile,line);
    if ( ((pos=line.find(olkey::olelse)) != std::string::npos) && (level==1) ) 
      trueclause=false;
    else if ( (pos=line.find(olkey::olendif)) != std::string::npos) 
      level--;
    else if ( !(trueclause ^ condition) ) // xor bitwise operator
      parse_oneline(line,infile);
    else {
      if ( (pos=line.find(olkey::iftrue)) != std::string::npos) 
	level++; // check for opening iftrue statement
      else if ( (pos=line.find(olkey::ifntrue)) != std::string::npos) 
	level++; // check for opening ifntrue statement
    }
  }
  return level?false:true ;
} 

void localSolverClient::parse_onefile(std::string fileName) { 
  int pos;
  std::string fullName=getWorkingDir()+fileName;
  std::ifstream infile(fullName.c_str());
  if (infile.is_open()){
    Msg::Info("Parse file <%s>",fullName.c_str());
    while (infile.good()){
      std::string line;
      getline(infile,line);
      parse_oneline(line,infile);
    }
    infile.close();
  }
  else
    Msg::Info("The file %s does not exist",fullName.c_str());
} 

bool localSolverClient::convert_ifstatement(std::ifstream &infile, std::ofstream &outfile, bool condition) { 
  int level, pos;
  std::string line;

  bool trueclause=true; 
  level=1;
  while ( infile.good() && level) {
    getline (infile,line);
    if ( ((pos=line.find(olkey::olelse)) != std::string::npos) && (level==1) ) 
      trueclause=false;
    else if ( (pos=line.find(olkey::olendif)) != std::string::npos) 
     level--;
    else if ( !(trueclause ^ condition) ) // xor bitwise operator
      convert_oneline(line,infile,outfile);
    else {
      if ( (pos=line.find(olkey::iftrue)) != std::string::npos) 
	level++; // check for opening iftrue statement
      else if ( (pos=line.find(olkey::ifntrue)) != std::string::npos) 
	level++; // check for opening ifntrue statement
    }
  }
  return level?false:true ;
} 

void localSolverClient::convert_oneline(std::string line, std::ifstream &infile, std::ofstream &outfile) { 
  int pos,cursor;
  std::vector<std::string> arguments;
  std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;
  std::vector<onelab::region> regions;
  std::string buff;

  if((pos=line.find_first_not_of(" \t"))==std::string::npos){
    // empty line
    outfile << line << std::endl;  
  }
  else if(!line.compare(pos,olkey::comment.size(),olkey::comment)){
    // commented out, skip the line
  }
  else if ( (pos=line.find(olkey::deflabel)) != std::string::npos){
    // onelab.tags(label,comment,separator)
    cursor = pos+olkey::deflabel.length();
    pos=line.find_first_of(')',cursor)+1;
    int NumArg=enclosed(line.substr(cursor,pos-cursor),arguments);
    if (NumArg == 0)
      Msg::Error("No onelab tags given");
    if(NumArg >= 1){
      olkey::label.assign(arguments[0]);
      olkey::line.assign(olkey::label+"line");
      olkey::begin.assign(olkey::label+"begin");
      olkey::end.assign(olkey::label+"end");
      olkey::include.assign(olkey::label+"include");
      olkey::ifcond.assign(olkey::label+"if");
      olkey::ifequal.assign(olkey::label+"ifequal");
      olkey::iftrue.assign(olkey::label+"iftrue");
      olkey::ifntrue.assign(olkey::label+"ifntrue"); 
      olkey::olelse.assign(olkey::label+"else");
      olkey::olendif.assign(olkey::label+"endif"); 
      olkey::getValue.assign(olkey::label+"getValue");
      olkey::getRegion.assign(olkey::label+"getRegion");
    }
    if(NumArg >= 2)
      olkey::comment.assign(arguments[1]);
    if(NumArg >= 3)
      olkey::separator.assign(arguments[2]);
    Msg::Info("Using now onelab tags <%s,%s,%s>",
	      olkey::label.c_str(),olkey::comment.c_str(),
	      olkey::separator.c_str());
  }
  else if( (pos=line.find(olkey::begin)) != std::string::npos) {
    // onelab.begin
    while (infile.good()){
      getline (infile,line);
      if( (pos=line.find(olkey::end)) != std::string::npos) return;
    }
    Msg::Fatal("Misformed <%s> block <%s>",
	       olkey::begin.c_str(),olkey::end.c_str());
  }
  else if ( (pos=line.find(olkey::iftrue)) != std::string::npos) {
    // onelab.iftrue
    cursor = pos+olkey::iftrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::iftrue.c_str(),line.c_str());
    bool condition = false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= true;
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      else
	Msg::Fatal("Unknown parameter <%s> in <%s> statement",
		   arguments[0].c_str(),olkey::iftrue.c_str());
    }
    if (!convert_ifstatement(infile,outfile,condition))
      Msg::Fatal("Misformed <%s> statement: %s",
		 olkey::iftrue.c_str(),arguments[0].c_str());     
  }
  else if ( (pos=line.find(olkey::ifntrue)) != std::string::npos) {
    // onelab.ifntrue
    cursor = pos+olkey::ifntrue.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::ifntrue.c_str(),line.c_str());
    bool condition = false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition= true;
    else{
      get(numbers,longName(arguments[0]));
      if (numbers.size())
	condition = (bool) numbers[0].getValue();
      else
	Msg::Fatal("Unknown parameter <%s> in <%s> statement",
		   arguments[0].c_str(),olkey::ifntrue.c_str());
    }
    if (!convert_ifstatement(infile,outfile,!condition))
      Msg::Fatal("Misformed <%s> statement: %s",
		 olkey::ifntrue.c_str(),arguments[0].c_str());  
  }
  else if ( (pos=line.find(olkey::ifcond)) != std::string::npos) {
    // onelab.ifcond
    cursor = pos+olkey::ifcond.length();
    int NumArgs=extractLogic(line.substr(cursor),arguments);
    bool condition= resolveLogicExpr(arguments);
    if (!convert_ifstatement(infile,outfile,condition))
      Msg::Fatal("Misformed %s statement: <%s>", line.c_str());
  }
  else if ( (pos=line.find(olkey::ifequal)) != std::string::npos) {
    // onelab.ifequal
    cursor = pos+olkey::ifequal.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<2)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::ifequal.c_str(),line.c_str());;
    bool condition= false;
    get(strings,longName(arguments[0]));
    if (strings.size())
      condition =  !strings[0].getValue().compare(arguments[1]);
    if (!convert_ifstatement(infile,outfile,condition))
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::ifequal.c_str(),line.c_str());
  }
  else if ( (pos=line.find(olkey::include)) != std::string::npos) {
    // onelab.include
    cursor = pos+olkey::include.length();
    pos=line.find_first_of(')',cursor)+1;
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::include.c_str(),line.c_str());
    convert_onefile(arguments[0],outfile);
  }
  else if ( (pos=line.find(olkey::getValue)) != std::string::npos) {
    // onelab.getValue
    // onelab.getValue, possibly several times in the line
    cursor=0;
    while ( (pos=line.find(olkey::getValue,cursor)) != std::string::npos){
      int pos0=pos; // for further use
      cursor = pos+olkey::getValue.length();
      pos=line.find_first_of(')',cursor)+1;
      if(enclosed(line.substr(cursor,pos-cursor),arguments)<1)
	Msg::Fatal("Misformed <%s> statement: (%s)",
		   olkey::getValue.c_str(),line.c_str());
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
  else if ( (pos=line.find(olkey::getRegion)) != std::string::npos) {
    // onelab.region
    int pos0=pos;
    cursor = pos+olkey::getRegion.length();
    pos=line.find_first_of(')',cursor)+1;
    int NumArg=enclosed(line.substr(cursor,pos-cursor),arguments);
    if(enclosed(line.substr(cursor,pos-cursor),arguments)<4)
      Msg::Fatal("Misformed <%s> statement: (%s)",
		 olkey::getRegion.c_str(),line.c_str());
    std::string paramName;
    paramName.assign(longName(arguments[0]));
    buff.assign(arguments[1]);
    std::cout << "getRegion:<" << arguments[0] << "> => <" << paramName << ">" << std::endl;
    get(regions,paramName);
    if (regions.size()){
      std::set<std::string> region;
      region=regions[0].getValue();
      for(std::set<std::string>::const_iterator it = region.begin();
          it != region.end(); it++){
	if(it != region.begin())  
	  buff.append(arguments[2].compare("comma")?arguments[2]:",");
        buff.append((*it));
      }
      buff.append(arguments[3]);
    }
    else{
      Msg::Fatal("Unknown region: <%s>",paramName.c_str());
    }
    line.replace(pos0,pos-pos0,buff); 
    outfile << line << std::endl; 
  }
  else if ( (pos=line.find(olkey::label)) != std::string::npos) {
    Msg::Warning("Ambiguous sentence: %s",line.c_str());
    Msg::Info("Using now onelab tags <%s,%s,%s>",olkey::label.c_str(),
	      olkey::comment.c_str(), olkey::separator.c_str());
  }
  else{
    outfile << line << std::endl; 
  }
}

void localSolverClient::convert_onefile(std::string fileName, std::ofstream &outfile) {
  int pos;
  std::string fullName=getWorkingDir()+fileName;
  std::ifstream infile(fullName.c_str());
  if (infile.is_open()){
    while ( infile.good() ) {
      std::string line;
      getline (infile,line);
      convert_oneline(line,infile,outfile);
    }
    infile.close();
  }
  else
    Msg::Fatal("The file %s cannot be opened",fullName.c_str());
}


void localSolverClient::client_sentence(const std::string &name, const std::string &action, 
		       const std::vector<std::string> &arguments) {
  Msg::Fatal("The action <%s> is unknown in this ccontext",action.c_str());
}

void MetaModel::client_sentence(const std::string &name, const std::string &action, 
		 const std::vector<std::string> &arguments){
  //std::vector<onelab::number> numbers;
  std::vector<onelab::string> strings;

  if(!action.compare("Register")){
    // syntax name.Register([interf...|encaps...]{,cmdl{,wdir,{host{,rdir}}}}) ;
    if(!findClientByName(name)){
      Msg::Info("Define client <%s>", name.c_str());

      std::string cmdl="",wdir="",host="",rdir="";
      if(arguments.size()>=2) cmdl.assign(resolveGetVal(arguments[1]));
      if(arguments.size()>=3) wdir.assign(resolveGetVal(arguments[2]));
      if(arguments.size()>=4) host.assign(resolveGetVal(arguments[3]));
      if(arguments.size()>=5) rdir.assign(resolveGetVal(arguments[4]));

      // check if one has a saved command line on the server 
      // (prealably read from file .ol.save)
      if(cmdl.empty())
	cmdl=Msg::GetOnelabString(name + "/CommandLine");
      registerClient(name,resolveGetVal(arguments[0]),cmdl,
		     getWorkingDir()+wdir,host,rdir);
      onelab::string o(name + "/CommandLine","");
      o.setValue(cmdl);
      o.setKind("file");
      o.setVisible(cmdl.empty());
      o.setAttribute("Highlight","Ivory");
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
	  c->setActive(atof( resolveGetVal(arguments[0]).c_str() ));
	else
	  Msg::Error("No argument for <%s.Active> statement", name.c_str());
      }
    }
    else
      Msg::Fatal("Unknown client <%s>", name.c_str());
  }
  else if(!action.compare(olkey::arguments)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/Arguments");
      strings[0].setValue(resolveGetVal(arguments[0]));
      strings[0].setVisible(false);
      set(strings[0]);
    }
  }
  else if(!action.compare(olkey::inFiles)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/InputFiles");
      strings[0].setValue(resolveGetVal(arguments[0]));
      strings[0].setKind("file");
      strings[0].setVisible(false);
      std::vector<std::string> choices;
      for(unsigned int i = 0; i < arguments.size(); i++)
	//if(std::find(choices.begin(),choices.end(),arguments[i])
	//==choices.end())
	choices.push_back(resolveGetVal(arguments[i]));
      strings[0].setChoices(choices);
      set(strings[0]);
    }
  }
  else if(!action.compare(olkey::outFiles)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/OutputFiles");
      strings[0].setValue(resolveGetVal(arguments[0]));
      strings[0].setKind("file");
      strings[0].setVisible(false);
      std::vector<std::string> choices;
      for(unsigned int i = 0; i < arguments.size(); i++)
	//if(std::find(choices.begin(),choices.end(),arguments[i])
	//==choices.end())
	choices.push_back(resolveGetVal(arguments[i]));
      strings[0].setChoices(choices);
      set(strings[0]);
    }
  }
  else if(!action.compare(olkey::upload)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/PostArray");
      strings[0].setValue(resolveGetVal(arguments[0]));
      std::vector<std::string> choices;
      for(unsigned int i = 0; i < arguments.size(); i++)
	choices.push_back(resolveGetVal(arguments[i]));
      strings[0].setChoices(choices);
      strings[0].setVisible(false);
      set(strings[0]);
    }
  }
  else if(!action.compare(olkey::merge)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/Merge");
      strings[0].setValue(resolveGetVal(arguments[0]));
      strings[0].setKind("file");
      strings[0].setVisible(false);
      std::vector<std::string> choices;
      for(unsigned int i = 0; i < arguments.size(); i++)
	choices.push_back(resolveGetVal(arguments[i]));
      strings[0].setChoices(choices);
      set(strings[0]);
    }
  }
  else if(!action.compare(olkey::checkCmd)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/9CheckCommand");
      strings[0].setValue(resolveGetVal(arguments[0]));
      strings[0].setVisible(false);
      set(strings[0]);
    }
  }
  else if(!action.compare(olkey::computeCmd)){
    if(arguments[0].size()){
      strings.resize(1);
      strings[0].setName(name+"/9ComputeCommand");
      strings[0].setValue(resolveGetVal(arguments[0]));
      strings[0].setVisible(false);
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
	if(std::find(choices.begin(),choices.end(),arguments[i])
	   ==choices.end())
	  choices.push_back(resolveGetVal(arguments[i]));
      strings[0].setChoices(choices);
      set(strings[0]);
    }
  }
  else if(!action.compare("List")){
    //no check whether choices[i] already inserted
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
    Msg::Fatal("Unknown action <%s>",action.c_str());
}
