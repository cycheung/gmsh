// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include "OnelabMessage.h"
#include "GmshSocket.h"
#include "OS.h"

#define ALWAYS_TRUE 1

int Msg::_commRank = 0;
int Msg::_commSize = 1;
int Msg::_verbosity = 4;
int Msg::_progressMeterStep = 10;
int Msg::_progressMeterCurrent = 0;
std::map<std::string, double> Msg::_timers;
int Msg::_warningCount = 0;
int Msg::_errorCount = 0;
GmshMessage *Msg::_callback = 0;
std::string Msg::_commandLine;
std::string Msg::_launchDate;
GmshClient *Msg::_client = 0;
onelab::client *Msg::_onelabClient = 0;
bool Msg::hasGmsh=false;

#if defined(HAVE_NO_VSNPRINTF)
static int vsnprintf(char *str, size_t size, const char *fmt, va_list ap)
{
  if(strlen(fmt) > size - 1){ // just copy the format
    strncpy(str, fmt, size - 1);
    str[size - 1] = '\0';
    return size;
  }
  return vsprintf(str, fmt, ap);
}
#endif

#if defined(_MSC_VER) && (_MSC_VER == 1310) //NET 2003
#define vsnprintf _vsnprintf
#endif

void Msg::Init(int argc, char **argv)
{
  time_t now;
  time(&now);
  _launchDate = ctime(&now);
  _launchDate.resize(_launchDate.size() - 1);
  _commandLine.clear();
  for(int i = 0; i < argc; i++){
    if(i) _commandLine += " ";
    _commandLine += argv[i];
  }
}

void Msg::Exit(int level)
{
  // delete the temp file
  //if(!_commRank) UnlinkFile(CTX::instance()->homeDir + CTX::instance()->tmpFileName);

  // exit directly on abnormal program termination (level != 0). We
  // used to call abort() to flush open streams, but on modern OSes
  // this calls the annoying "report this crash to the mothership"
  // window... so just exit!
  if(level){
    exit(level);
  }
  exit(_errorCount);
}

void Msg::Fatal(const char *fmt, ...)
{
  _errorCount++;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback) (*_callback)("Fatal", str);
  if(_client) _client->Error(str);

  if(ALWAYS_TRUE){
    if(_commSize > 1)
      fprintf(stderr, "Fatal   : [On processor %d] %s\n", _commRank, str);
    else
      fprintf(stderr, "Fatal   : %s\n", str);
   fflush(stderr);
  }

  FinalizeClient();
  FinalizeOnelab();
  delete loader;
  // only exit if a callback is not provided
  //if(!_callback) Exit(1);
  Exit(1);
}

void Msg::Error(const char *fmt, ...)
{
  _errorCount++;

  if(_verbosity < 1) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback) (*_callback)("Error", str);
  if(_client) _client->Error(str);

  if(ALWAYS_TRUE){
    if(_commSize > 1)
      fprintf(stderr, "Error   : [On processor %d] %s\n", _commRank, str);
    else
      fprintf(stderr, "Error   : %s\n", str);
    fflush(stderr);
  }
}

void Msg::Warning(const char *fmt, ...)
{
  _warningCount++;

  if(_commRank || _verbosity < 2) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback) (*_callback)("Warning", str);
  if(_client) _client->Warning(str);

  if(ALWAYS_TRUE){
    fprintf(stderr, "Warning : %s\n", str);
    fflush(stderr);
  }
}

void Msg::Info(const char *fmt, ...)
{
  if(_commRank || _verbosity < 3) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback) (*_callback)("Info", str);
  if(_client) _client->Info(str);

  if(ALWAYS_TRUE){
    fprintf(stdout, "Info    : %s\n", str);
    fflush(stdout);
  }
}

void Msg::Direct(const char *fmt, ...)
{
  if(_commRank || _verbosity < 3) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  Direct(3, str);
}

void Msg::Direct(int level, const char *fmt, ...)
{
  if(_commRank || _verbosity < level) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback) (*_callback)("Direct", str);
  if(_client) _client->Info(str);

  if(ALWAYS_TRUE){
    fprintf(stdout, "%s\n", str);
    fflush(stdout);
  }
}

void Msg::StatusBar(int num, bool log, const char *fmt, ...)
{
  if(_commRank || _verbosity < 3) return;
  if(num < 1 || num > 3) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback && log) (*_callback)("Info", str);
  if(_client && log) _client->Info(str);

  if(log && ALWAYS_TRUE){
    fprintf(stdout, "Info    : %s\n", str);
    fflush(stdout);
  }
}

void Msg::Debug(const char *fmt, ...)
{
  if(_verbosity < 99) return;

  char str[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(str, sizeof(str), fmt, args);
  va_end(args);

  if(_callback) (*_callback)("Debug", str);
  if(_client) _client->Info(str);

  if(ALWAYS_TRUE){
    if(_commSize > 1)
      fprintf(stdout, "Debug   : [On processor %d] %s\n", _commRank, str);
    else
      fprintf(stdout, "Debug   : %s\n", str);
    fflush(stdout);
  }
}

void Msg::ProgressMeter(int n, int N, const char *fmt, ...)
{
  if(_commRank || _verbosity < 3) return;

  double percent = 100. * (double)n/(double)N;

  if(percent >= _progressMeterCurrent){
    char str[1024];
    va_list args;
    va_start(args, fmt);
    vsnprintf(str, sizeof(str), fmt, args);
    va_end(args);

    if(strlen(fmt)) strcat(str, " ");

    char str2[1024];
    sprintf(str2, "(%d %%)", _progressMeterCurrent);
    strcat(str, str2);

    if(_client) _client->Progress(str);

    if(ALWAYS_TRUE){
      fprintf(stdout, "%s                     \r", str);
      fflush(stdout);
    }

    while(_progressMeterCurrent < percent)
      _progressMeterCurrent += _progressMeterStep;
  }

  if(n > N - 1){
    if(_client) _client->Progress("Done!");

    if(ALWAYS_TRUE){
      fprintf(stdout, "Done!                                              \r");
      fflush(stdout);
    }
  }
}

void Msg::PrintTimers()
{
  // do a single stdio call!
  std::string str;
  for(std::map<std::string, double>::iterator it = _timers.begin();
      it != _timers.end(); it++){
    if(it != _timers.begin()) str += ", ";
    char tmp[256];
    sprintf(tmp, "%s = %gs ", it->first.c_str(), it->second);
    str += std::string(tmp);
  }
  if(!str.size()) return;

  if(ALWAYS_TRUE){
    if(_commSize > 1)
      fprintf(stdout, "Timers  : [On processor %d] %s\n", _commRank, str.c_str());
    else
      fprintf(stdout, "Timers  : %s\n", str.c_str());
    fflush(stdout);
  }
}

void Msg::PrintErrorCounter(const char *title)
{
  if(_commRank || _verbosity < 1) return;
  if(!_warningCount && !_errorCount) return;

  std::string prefix = _errorCount ? "Error   : " : "Warning : ";
  std::string help("Check the full log for details");
  std::string line(std::max(strlen(title), help.size()), '-');
  char warn[128], err[128];
  sprintf(warn, "%5d warning%s", _warningCount, _warningCount == 1 ? "" : "s");
  sprintf(err, "%5d error%s", _errorCount, _errorCount == 1 ? "" : "s");

  if(ALWAYS_TRUE){
    fprintf(stderr, "%s\n%s\n%s\n%s\n%s\n%s\n", (prefix + line).c_str(),
            (prefix + title).c_str(), (prefix + warn).c_str(),
            (prefix + err).c_str(), (prefix + help).c_str(),
            (prefix + line).c_str());
    fflush(stderr);
  }
}

double Msg::GetValue(const char *text, double defaultval)
{
  printf("%s (default=%.16g): ", text, defaultval);
  char str[256];
  char *ret = fgets(str, sizeof(str), stdin);
  if(!ret || !strlen(str) || !strcmp(str, "\n"))
    return defaultval;
  else
    return atof(str);
}

std::string Msg::GetString(const char *text, std::string defaultval)
{
  printf("%s (default=%s): ", text, defaultval.c_str());
  char str[256];
  char *ret = fgets(str, sizeof(str), stdin);
  if(!ret || !strlen(str) || !strcmp(str, "\n"))
    return defaultval;
  else
    return std::string(str);
}

int Msg::GetAnswer(const char *question, int defaultval, const char *zero,
                   const char *one, const char *two)
{
  if(two)
    printf("%s\n\n0=[%s] 1=[%s] 2=[%s] (default=%d): ", question,
           zero, one, two, defaultval);
  else
    printf("%s\n\n0=[%s] 1=[%s] (default=%d): ", question,
           zero, one, defaultval);
  char str[256];
  char *ret = fgets(str, sizeof(str), stdin);
  if(!ret || !strlen(str) || !strcmp(str, "\n"))
    return defaultval;
  else
    return atoi(ret);
}

void Msg::InitClient(std::string sockname)
{
  if(_client) delete _client;
  _client = new GmshClient();
  if(_client->Connect(sockname.c_str()) < 0){
    Msg::Error("Unable to connect to server on %s", sockname.c_str());
    delete _client;
    _client = 0;
  }
  else
    _client->Start();
}

void Msg::FinalizeClient()
{
  if(_client){
    _client->Stop();
    _client->Disconnect();
    delete _client;
  }
  _client = 0;
}

void Msg::Barrier()
{
}

void Msg::InitializeOnelab(const std::string &name, const std::string &sockname)
{
  if(_onelabClient) delete _onelabClient;
  if (sockname.empty())
    _onelabClient = new onelab::localClient(name);
  else{
    onelab::remoteNetworkClient *c = new onelab::remoteNetworkClient(name, sockname);
    _onelabClient = c;
    _client = c->getGmshClient();
  }
}

double Msg::GetOnelabNumber(std::string name)
{
  if(_onelabClient){
    std::vector<onelab::number> ps;
    _onelabClient->get(ps, name);
    if(ps.size())
      return ps[0].getValue();
  }
  return 0;
}

void Msg::GetOnelabNumber(std::string name, double *val)
{
  if(_onelabClient){
    std::vector<onelab::number> ps;
    _onelabClient->get(ps, name);
    if(ps.size()){
      *val = ps[0].getValue();
      return;
    }
  }
  *val = 0;
}

void Msg::SetOnelabNumber(std::string name, double val, bool visible)
{
  if(_onelabClient){
    onelab::number o(name, val);
    o.setVisible(visible);
    _onelabClient->set(o);
  }
}

void Msg::SetOnelabNumber(onelab::number s)
{
  if(_onelabClient){
    _onelabClient->set(s);
  }
}

std::string Msg::GetOnelabString(std::string name)
{
  std::string str="";
  if(_onelabClient){
    std::vector<onelab::string> ps;
    _onelabClient->get(ps, name);
    if(ps.size() && ps[0].getValue().size())
      str = ps[0].getValue();
  }
  return str;
}

bool Msg::GetOnelabChoices(std::string name, std::vector<std::string> &choices){
  if(_onelabClient){
    std::vector<onelab::string> ps;
    _onelabClient->get(ps, name);
    if(ps.size() && ps[0].getValue().size()){
      choices=ps[0].getChoices();
      return true;
    }
  }
  return false;
}

void Msg::SetOnelabString(std::string name, std::string val, bool visible)
{
  if(_onelabClient){
    onelab::string o(name, val);
    o.setVisible(visible);
    _onelabClient->set(o);
  }
  else
    std::cout << "Pas de client" << std::endl;
}

void Msg::SetOnelabString(onelab::string s){
  if(_onelabClient){
    _onelabClient->set(s);
  }
}

void Msg::SetOnelabRegion(onelab::region r){
  if(_onelabClient){
    _onelabClient->set(r);
  }
}

void Msg::SetOnelabAttributeString(std::string name,
				   std::string attrib,std::string val){
  if(_onelabClient){
    std::vector<onelab::string> ps;
    _onelabClient->get(ps, name);
    if(ps.size()){
      ps[0].setAttribute(attrib,val);
    }
  }
}
std::string Msg::GetOnelabAttributeString(std::string name,std::string attrib){
  std::string str="";
  if(_onelabClient){
    std::vector<onelab::string> ps;
    _onelabClient->get(ps, name);
    if(ps.size())
      str = ps[0].getAttribute(attrib);
  }
  return str;
}
std::string Msg::GetOnelabAttributeNumber(std::string name,std::string attrib){
  std::string str="";
  if(_onelabClient){
    std::vector<onelab::number> ps;
    _onelabClient->get(ps, name);
    if(ps.size())
      str = ps[0].getAttribute(attrib);
  }
  return str;
}

/* not used 
void Msg::ExchangeOnelabParameter(const std::string &key,
                                  std::vector<double> &val,
                                  std::map<std::string, std::vector<double> > &fopt,
                                  std::map<std::string, std::vector<std::string> > &copt)
{
  if(!_onelabClient || val.empty()) return;

  std::string name(key);
  if(copt.count("Path")){
    std::string path = copt["Path"][0];
    // if path ends with a number, assume it's for ordering purposes
    if(path.size() && path[path.size() - 1] >= '0' && path[path.size() - 1] <= '9')
      name = path + name;
    else if(path.size() && path[path.size() - 1] == '/')
      name = path + name;
    else
      name = path + "/" + name;
  }
  std::vector<onelab::number> ps;
  _onelabClient->get(ps, name);
  if(ps.size()){ // use value from server
    val[0] = ps[0].getValue();
  }
  else{ // send value to server
    onelab::number o(name, val[0]);
    if(fopt.count("Range") && fopt["Range"].size() == 2){
      o.setMin(fopt["Range"][0]); o.setMax(fopt["Range"][1]);
    }
    else if(fopt.count("Min") && fopt.count("Max")){
      o.setMin(fopt["Min"][0]); o.setMax(fopt["Max"][0]);
    }
    else if(fopt.count("Min")){
      o.setMin(fopt["Min"][0]); o.setMax(1.e200);
    }
    else if(fopt.count("Max")){
      o.setMax(fopt["Max"][0]); o.setMin(-1.e200);
    }
    if(fopt.count("Step")) o.setStep(fopt["Step"][0]);
    if(fopt.count("Choices")) o.setChoices(fopt["Choices"]);
    if(copt.count("Help")) o.setHelp(copt["Help"][0]);
    if(copt.count("Label")) o.setLabel(copt["Label"][0]);
    _onelabClient->set(o);
  }
}
*/

void Msg::AddOnelabNumberChoice(std::string name, double val)
{
  if(_onelabClient){
    std::vector<double> choices;
    std::vector<onelab::number> ps;
    _onelabClient->get(ps, name);
    if(ps.size()){
      choices = ps[0].getChoices();
    }
    else{
      ps.resize(1);
      ps[0].setName(name);
    }
    ps[0].setValue(val);
    choices.push_back(val);
    ps[0].setChoices(choices);
    ps[0].setAttribute("Highlight","Coral"); // only used by PostArray
    ps[0].setReadOnly(true);
    _onelabClient->set(ps[0]);
  }
}

void Msg::AddOnelabStringChoice(std::string name, std::string kind,
                                    std::string value)
{
  if(_onelabClient){
    std::vector<std::string> choices;
    std::vector<onelab::string> ps;
    _onelabClient->get(ps, name);
    if(ps.size()){
      choices = ps[0].getChoices();
      if(std::find(choices.begin(), choices.end(), value) == choices.end())
        choices.push_back(value);
    }
    else{
      ps.resize(1);
      ps[0].setName(name);
      ps[0].setKind(kind);
      choices.push_back(value);
    }
    ps[0].setValue(value);
    ps[0].setChoices(choices);
    _onelabClient->set(ps[0]);
  }
}

int Msg::Synchronize_Down(){
  std::vector<onelab::number> numbers;
  onelab::number *pn;
  loader->get(numbers,"");
  if(numbers.size()){
    for(std::vector<onelab::number>::const_iterator it = numbers.begin();
	it != numbers.end(); it++){
      pn=new onelab::number;
      pn->fromChar((*it).toChar());
      Msg::SetOnelabNumber(*pn);
      delete pn;
    }
  }
  std::vector<onelab::string> strings;
  onelab::string *ps;
  loader->get(strings,"");
  if(strings.size()){
    for(std::vector<onelab::string>::const_iterator it = strings.begin();
  	it != strings.end(); it++){
      ps=new onelab::string;
      ps->fromChar((*it).toChar());
      Msg::SetOnelabString(*ps);
      delete ps;
    }
  }
  std::vector<onelab::region> regions;
  onelab::region *pr;
  loader->get(regions,"");
  if(regions.size()){
    for(std::vector<onelab::region>::const_iterator it = regions.begin();
  	it != regions.end(); it++){
      pr=new onelab::region;
      pr->fromChar((*it).toChar());
      Msg::SetOnelabRegion(*pr);
      delete pr;
    }
  }
  return(numbers.size()+strings.size()+regions.size());
}

int Msg::Synchronize_Up(){
  std::vector<onelab::number> numbers;
  onelab::number *pn;
  _onelabClient->get(numbers,"");
  if(numbers.size()){
    for(std::vector<onelab::number>::const_iterator it = numbers.begin();
  	it != numbers.end(); it++){
      pn = new(onelab::number);
      pn->fromChar((*it).toChar());
      loader->set(*pn);
      delete pn;
    }
  }
  std::vector<onelab::string> strings;
  onelab::string *ps;
  _onelabClient->get(strings,"");
  if(strings.size()){
    for(std::vector<onelab::string>::const_iterator it = strings.begin();
	it != strings.end(); it++){
      ps=new onelab::string;
      ps->fromChar((*it).toChar());
      loader->set(*ps);
      delete ps;
    }
  }
  std::vector<onelab::region> regions;
  onelab::region *pr;
  _onelabClient->get(regions,"");
  if(regions.size()){
    for(std::vector<onelab::region>::const_iterator it = regions.begin();
	it != regions.end(); it++){
      pr=new onelab::region;
      pr->fromChar((*it).toChar());
      loader->set(*pr);
      delete pr;
    }
  }
  return(numbers.size()+strings.size()+regions.size());
}

void Msg::FinalizeOnelab(){
  if(_onelabClient){
    delete _onelabClient;
    _onelabClient = 0;
    _client = 0;
  }
}


