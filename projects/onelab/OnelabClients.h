// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _ONELAB_CLIENTS_H_
#define _ONELAB_CLIENTS_H_

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include "OS.h"
#include "onelab.h"
#include "OnelabMessage.h"

namespace olkey{ // reserved keywords for onelab
  static std::string label("OL.");
  static std::string client(label+"client");
  static std::string param(label+"parameter");
  static std::string number(label+"number"), string(label+"string");
  static std::string include(label+"include"); 
  static std::string iftrue(label+"iftrue"), olelse(label+"else"), olendif(label+"endif"); 
  static std::string ifequal(label+"ifequal");
  static std::string getValue(label+"getValue");
  static std::string extension(".ol");
  static std::string comment("%");
  static std::string arguments("Args"), inFiles("In"), outFiles("Out");
}

static char charSep() { return '\0'; }
#if defined(WIN32)
//static char[1] dirSep() { return "\"; }
static std::string dirSep("\\");
#else
static std::string dirSep("/");
//static char[1] dirSep() { return "/"; }
#endif


int getOptions(int argc, char *argv[], std::string &action, std::string &commandLine, std::string &fileName, std::string &clientName, std::string &sockName, int &modelNumber);
std::string itoa(const int i);
//bool fileExist(std::string filename);
//int newStep();
void appendOption(std::string &str, const std::string &what, const int val);
void appendOption(std::string &str, const std::string &what);
void GmshDisplay(onelab::remoteNetworkClient *loader, std::string fileName, std::vector<std::string> choices);
//void GmshDisplay(onelab::remoteNetworkClient *loader, std::string modelName, std::string fileName);
std::string getCurrentWorkdir();
std::string getUserHomedir();
static std::string getNextToken(const std::string &msg,std::string::size_type &first);
std::string sanitize(const std::string &in);
std::string removeBlanks(const std::string &in);
int enclosed(const std::string &in, std::vector<std::string> &arguments);
int extract(const std::string &in, std::string &paramName, std::string &action, std::vector<std::string> &arguments);
bool extractRange(const std::string &in, std::vector<double> &arguments);

bool checkIfPresent(std::string filename);

typedef std::vector <std::vector <double> > array;
array read_array(std::string filename, char sep);
double find_in_array(const int i, const int j, const std::vector <std::vector <double> > &data);
std::vector<double> extract_column(const int j, array data);

class PromptUser : public onelab::localClient {
public:
 PromptUser(const std::string &name) : onelab::localClient(name) {}
  ~PromptUser(){}
  int  getVerbosity();
  void setVerbosity(const int ival);
  void setNumber(const std::string paramName, const double val, const std::string &help="");
  double getNumber(const std::string paramName);
  bool existNumber(const std::string paramName);
  void setString(const std::string paramName, const std::string &val, const std::string &help="");
  std::string getString(const std::string paramName);
  bool existString(const std::string paramName);
  void addNumberChoice(std::string name, double val);
  void addStringChoice(std::string name, std::string str);
  //std::vector<std::string> getChoices(const std::string paramName);
  std::string stateToChar();  
  std::string showParamSpace();
  std::string showClientStatus();
  void statusClients();
  bool menu(std::string commandLine, std::string fileName, int modelNumber);
};

static std::string getShortName(const std::string &name) {
  std::string s = name;
  // remove path
  std::string::size_type last = name.find_last_of('/');
  if(last != std::string::npos)
    s = name.substr(last + 1);
  // remove starting numbers
  while(s.size() && s[0] >= '0' && s[0] <= '9')
    s = s.substr(1);
  return s;
}

class ShortNameLessThan{
 public:
  bool operator()(const std::string p1, const std::string p2) const
  {
    return getShortName(p1) < getShortName(p2);
  }
};

/*

VIRTUAL and BASE CLASSES

localSolverClient est la classe de base pour tous les clients
avec les méthodes (virtuelles) analyze() et compute()  (que la classe 'localClient' n'a pas)
qui sont les deux modes d'exécution du métamodèle.
Seule _commandLine est stockée dans la classe
Les autres infos sont définies sur le serveur
*/
class localSolverClient : public onelab::localClient{
 private:
  std::string _commandLine;
  bool _enabled;
  std::set<std::string, ShortNameLessThan> _parameters;
  std::string longName(const std::string name);

  std::string evaluateGetVal(std::string line);
 public:
 localSolverClient(const std::string &name, const std::string &commandLine) 
   : onelab::localClient(name), _commandLine(commandLine), _enabled(true) {
  }
  virtual ~localSolverClient(){}
  const std::string &getCommandLine(){ return _commandLine; }
  virtual void setCommandLine(const std::string &s){ _commandLine = s; }
  const std::string getString(const std::string what);
  const bool getList(const std::string type, std::vector<std::string> &choices);
  const bool getActive() { return _enabled; }
  const void setActive(int val) { _enabled=(bool)val; }
  const std::string buildRmCommand();

  virtual bool checkCommandLine();
  virtual std::string toChar();

  std::string resolveGetVal(std::string line);
  void parse_onefile(std::string ifilename);
  virtual void parse_clientline(std::string line, std::ifstream &infile) {}
  void parse_oneline(std::string line, std::ifstream &infile) ;
  bool parse_ifstatement(std::ifstream &infile, bool condition) ;
  void convert_onefile(std::string ifilename, std::ofstream &outfile);
  void convert_oneline(std::string line, std::ifstream &infile, std::ofstream &outfile);
  bool convert_ifstatement(std::ifstream &infile, std::ofstream &outfile, bool condition) ;

  virtual void analyze() =0;
  virtual void compute() =0;
};

class localNetworkSolverClient : public localSolverClient{
 private:
  // command line option to specify socket
  std::string _socketSwitch;
  // pid of the remote network client
  int _pid;
  // underlying GmshServer
  GmshServer *_gmshServer;
 public:
 localNetworkSolverClient(const std::string &name, const std::string &commandLine)
   : localSolverClient(name,commandLine), _socketSwitch("-onelab"),
    _pid(-1), _gmshServer(0) {}
  virtual ~localNetworkSolverClient(){}
  virtual bool isNetworkClient(){ return true; }
  const std::string &getSocketSwitch(){ return _socketSwitch; }
  void setSocketSwitch(const std::string &s){ _socketSwitch = s; }
  int getPid(){ return _pid; }
  void setPid(int pid){ _pid = pid; }
  GmshServer *getGmshServer(){ return _gmshServer; }
  void setGmshServer(GmshServer *server){ _gmshServer = server; }

  virtual std::string buildCommandLine();
  virtual bool run();
  virtual bool kill();

  virtual void analyze() =0;
  virtual void compute() =0;
};

class remoteClient {
 private:
  std::string _remoteHost, _remoteDir;
 public:
 remoteClient(const std::string &host, const std::string &dir) 
   : _remoteHost(host), _remoteDir(dir) {}
  ~remoteClient(){}

  const std::string &getRemoteHost() const { return _remoteHost; }
  const std::string &getRemoteDir() const { return _remoteDir; }

  bool checkIfPresentRemote(const std::string &fileName);
  bool syncInputFile(const std::string &fileName);
  bool syncOutputFile(const std::string &fileName);
};

// ONELAB CLIENTS

class MetaModel : public localSolverClient {
private:
  std::vector<localSolverClient *> _clients;
  void registerInterfacedClient(const std::string name, const std::string cmdl);
  void registerDistantClient(const std::string name, const std::string cmdl, const std::string host, const std::string dir);
  void registerEncapsulatedClient(const std::string name, const std::string cmdl);
 public:
 MetaModel(const std::string &commandLine, const std::string &cname, const std::string &fname, const int number) 
   : localSolverClient(cname,commandLine){
    clientName = cname;
    modelNumberFromArgs = number;
    genericNameFromArgs = fname.size() ? fname : commandLine;
    parse_onefile(genericNameFromArgs + olkey::extension);
    parse_onefile(genericNameFromArgs + olkey::extension + ".save");
  }
  ~MetaModel(){}
  typedef std::vector<localSolverClient*>::iterator citer;
  citer firstClient(){ return _clients.begin(); }
  citer lastClient(){ return _clients.end(); }
  int getNumClients() { return _clients.size(); };

  void registerClient(const std::string &name, const std::string &type, const std::string &cmdl);
  void registerClient(const std::string &name, const std::string &type, const std::string &cmdl, 
		      const std::string &host, const std::string &dir);

  bool checkCommandLines();
  void saveCommandLines(const std::string fileName);
  localSolverClient *findClientByName(std::string name){
    for(unsigned int i=0; i<_clients.size(); i++)
      if(_clients[i]->getName() == name) return _clients[i];
    return 0;
  }
  std::string genericNameFromArgs, clientName;
  int modelNumberFromArgs;
  void parse_clientline(std::string line, std::ifstream &infile);
  std::string toChar(){}
  void PostArray(std::vector<std::string> choices);
  void initialize();
  void simpleCheck();
  void simpleCompute();
  void analyze();
  void compute();
};

class InterfacedClient : public localSolverClient { 
  // n'utilise pas localNetworkSolverClient::run mais client::run()
 public:
 InterfacedClient(const std::string &name, const std::string &commandLine)
   : localSolverClient(name,commandLine) {}
  ~InterfacedClient(){}

  void analyze();
  void convert();
  virtual void compute();
};

class EncapsulatedClient : public localNetworkSolverClient { 
  // utilise localNetworkClient::run
public:
 EncapsulatedClient(const std::string &name, const std::string &commandLine) 
   : localNetworkSolverClient(name,commandLine) {}
  ~EncapsulatedClient(){}

  void analyze();
  void compute() ;
};

class RemoteInterfacedClient : public InterfacedClient, public remoteClient {
public:
 RemoteInterfacedClient(const std::string &name, const std::string &commandLine, const std::string &host, const std::string &dir) 
   : InterfacedClient(name,commandLine), remoteClient(host,dir) {}
  ~RemoteInterfacedClient(){}

  bool checkCommandLine();
  /* void analyze(); */
  /* void convert(); */
  void compute() ;
};

class RemoteEncapsulatedClient : public EncapsulatedClient, public remoteClient {
public:
 RemoteEncapsulatedClient(const std::string &name, const std::string &commandLine, const std::string &host, const std::string &dir) 
   : EncapsulatedClient(name,commandLine), remoteClient(host,dir) {}
  ~RemoteEncapsulatedClient(){}

  std::string buildCommandLine();
  bool checkCommandLine();
  void compute() ;
};

#endif


