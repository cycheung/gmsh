// ONELAB - Copyright (C) 2010-2012 C. Geuzaine, F. Henrotte
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

// Onelab file extension
static std::string onelabExtension(".ol");

static char charSep() { return '\0'; }
#if defined(WIN32)
static std::string dirSep("\\");
static std::string cmdSep(" & ");
#else
static std::string dirSep("/");
static std::string cmdSep(" ; ");
#endif

int getOptions(int argc, char *argv[], std::string &action, std::string &commandLine, std::string &caseName, std::string &clientName, std::string &sockName, int &modelNumber);
std::string itoa(const int i);
std::string ftoa(const double x);
void appendOption(std::string &str, const std::string &what, const int val);
void appendOption(std::string &str, const std::string &what);
void GmshDisplay(onelab::remoteNetworkClient *loader, std::string fileName, std::vector<std::string> choices);
std::string getCurrentWorkdir();
std::string getUserHomedir();
//static std::string getNextToken(const std::string &msg,std::string::size_type &first);
std::string sanitize(const std::string &in);
std::string removeBlanks(const std::string &in);
int enclosed(const std::string &in, std::vector<std::string> &arguments);
int extract(const std::string &in, std::string &paramName, std::string &action, std::vector<std::string> &arguments);
bool extractRange(const std::string &in, std::vector<double> &arguments);
std::string extractExpandPattern(const std::string& str);
bool checkIfPresent(std::string filename);

typedef std::vector <std::vector <double> > array;
array read_array(std::string filename, char sep);
double find_in_array(int i, int j, const std::vector <std::vector <double> > &data);
std::vector<double> extract_column(const int j, array data);

class PromptUser : public onelab::localClient {
public:
 PromptUser(const std::string &name) : onelab::localClient(name) {}
  ~PromptUser(){}
  /* int  getVerbosity(); */
  /* void setVerbosity(const int ival); */
  void setNumber(const std::string paramName, const double val, const std::string &help="");
  double getNumber(const std::string paramName);
  bool existNumber(const std::string paramName);
  void setString(const std::string paramName, const std::string &val, const std::string &help="");
  std::string getString(const std::string paramName);
  bool existString(const std::string paramName);
  void addNumberChoice(std::string name, double val);
  void addStringChoice(std::string name, std::string str);
  std::string stateToChar();  
  std::string showParamSpace();
  std::string showClientStatus();
  void statusClients();
  bool menu(std::string commandLine, std::string workingDir, std::string fileName, int modelNumber);
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

localSolverClient est la classe de base pour tous les clients des métamodèles
Elle a les méthodes (virtuelles) analyze() et compute() 
(que la classe 'localClient' n'a pas)
qui sont les deux modes d'exécution du métamodèle.
Seule _commandLine et _workingDir sont stockées dans la classe
Les autres données décrivant le client 
(input et output files, arguments) sont stockées sur le serveur... 
*/
class localSolverClient : public onelab::localClient{
 private:
  std::string _commandLine;
  std::string _workingDir;
  int _active;
  bool _onelabBlock;
  std::set<std::string, ShortNameLessThan> _parameters;
  std::string longName(const std::string name);
  //std::string evaluateGetVal(std::string line);
 public:
 localSolverClient(const std::string &name, const std::string &cmdl, 
		   const std::string &wdir) 
   : onelab::localClient(name), _commandLine(cmdl), _workingDir(wdir),
    _active(1), _onelabBlock(false) {
  }
  virtual ~localSolverClient(){}
  const std::string &getCommandLine(){ return _commandLine; }
  const std::string &getWorkingDir() { return _workingDir; }
  virtual void setCommandLine(const std::string &s){ _commandLine = s; }
  virtual void setWorkingDir(const std::string &s){ _workingDir = s; }

  const std::string getString(const std::string what);
  const bool getList(const std::string type, 
		     std::vector<std::string> &choices);
  const bool isActive() { return (bool)_active; }
  const void setActive(int val) { _active=val; }
  int getActive() { return _active; }
  const bool isOnelabBlock() { return _onelabBlock; }
  const void openOnelabBlock() { _onelabBlock=true; }
  const void closeOnelabBlock() { _onelabBlock=false; }
  bool buildRmCommand(std::string &cmd);
  bool checkIfPresentLocal(const std::string &fileName){
    return checkIfPresent(getWorkingDir()+fileName);
  }
  virtual bool checkCommandLine();
  virtual std::string toChar();

  std::string resolveGetVal(std::string line);
  bool resolveLogicExpr(std::vector<std::string> arguments);
  virtual void client_sentence(const std::string &name, 
			       const std::string &action, 
			       const std::vector<std::string> &arguments);
  void modify_tags(const std::string lab, const std::string com);
  void parse_sentence(std::string line) ;
  void parse_oneline(std::string line, std::ifstream &infile) ;
  bool parse_block(std::ifstream &infile) ;
  bool parse_ifstatement(std::ifstream &infile, bool condition) ;
  void parse_onefile(std::string ifilename, bool mandatory=true);
  void convert_oneline(std::string line, std::ifstream &infile, 
		       std::ofstream &outfile);
  bool convert_ifstatement(std::ifstream &infile, 
			   std::ofstream &outfile, bool condition) ;
  void convert_onefile(std::string ifilename, std::ofstream &outfile);
  virtual void analyze() =0;
  virtual void compute() =0;
  void PostArray(std::vector<std::string> choices);
  void GmshMerge(std::vector<std::string> choices);
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
 localNetworkSolverClient(const std::string &name, const std::string &cmdl, const std::string &wdir)
   : localSolverClient(name,cmdl,wdir), _socketSwitch("-onelab"),
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
  std::string _remoteHost;
  std::string _remoteDir;
 public:
 remoteClient(const std::string &host, const std::string &rdir) 
   : _remoteHost(host), _remoteDir(rdir) {}
  ~remoteClient(){}

  const std::string &getRemoteHost() const { return _remoteHost; }
  const std::string &getRemoteDir() const { return _remoteDir; }

  bool checkIfPresentRemote(const std::string &fileName);
  bool syncInputFile(const std::string &wdir, const std::string &fileName);
  bool syncOutputFile(const std::string &wdir, const std::string &fileName);
};

// ONELAB CLIENTS

class MetaModel : public localSolverClient {
 private:
  std::vector<localSolverClient *> _clients;
 public:
 MetaModel(const std::string &cmdl, const std::string &wdir, const std::string &cname, const std::string &fname, const int number) 
   : localSolverClient(cname,cmdl,wdir){
    clientName = cname;
    modelNumberFromArgs = number;
    genericNameFromArgs = fname.size() ? fname : cmdl;
    openOnelabBlock();
    parse_onefile( genericNameFromArgs + onelabExtension + ".save",false);
    parse_onefile( genericNameFromArgs + onelabExtension);
    closeOnelabBlock();
  }
  ~MetaModel(){}
  typedef std::vector<localSolverClient*>::iterator citer;
  citer firstClient(){ return _clients.begin(); }
  citer lastClient(){ return _clients.end(); }
  int getNumClients() { return _clients.size(); };

  void registerClient(const std::string &name, const std::string &type, 
		      const std::string &cmdl, const std::string &wdir, 
		      const std::string &host, const std::string &rdir);
  bool checkCommandLines();
  void saveCommandLines(const std::string fileName);
  localSolverClient *findClientByName(std::string name){
    for(unsigned int i=0; i<_clients.size(); i++)
      if(_clients[i]->getName() == name) return _clients[i];
    return 0;
  }
  std::string genericNameFromArgs, clientName;
  int modelNumberFromArgs;
  void client_sentence(const std::string &name, const std::string &action, 
		       const std::vector<std::string> &arguments);
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
 InterfacedClient(const std::string &name, const std::string &cmdl, const std::string &wdir)
   : localSolverClient(name,cmdl,wdir) {}
  ~InterfacedClient(){}

  void analyze();
  void convert();
  virtual void compute();
};

class EncapsulatedClient : public localNetworkSolverClient { 
  // utilise localNetworkClient::run
public:
 EncapsulatedClient(const std::string &name, const std::string &cmdl, const std::string &wdir) 
   : localNetworkSolverClient(name,cmdl,wdir) {}
  ~EncapsulatedClient(){}

  void analyze();
  void compute() ;
};

class RemoteInterfacedClient : public InterfacedClient, public remoteClient {
public:
 RemoteInterfacedClient(const std::string &name, const std::string &cmdl, const std::string &wdir, const std::string &host, const std::string &rdir) 
   : InterfacedClient(name,cmdl,wdir), remoteClient(host,rdir) {}
  ~RemoteInterfacedClient(){}

  bool checkCommandLine();
  void compute() ;
};

class RemoteEncapsulatedClient : public EncapsulatedClient, public remoteClient {
public:
 RemoteEncapsulatedClient(const std::string &name, const std::string &cmdl, const std::string &wdir, const std::string &host, const std::string &rdir) 
   : EncapsulatedClient(name,cmdl,wdir), remoteClient(host,rdir) {}
  ~RemoteEncapsulatedClient(){}

  std::string buildCommandLine();
  bool checkCommandLine();
  void compute() ;
};

#endif


