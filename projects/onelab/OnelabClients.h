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

int getOptions(int argc, char *argv[], std::string &action, std::string &commandLine, std::string &fileName, std::string &clientName, std::string &sockName, int &modelNumber);
std::string itoa(const int i);
bool checkIfPresent(std::string filename);
bool checkIfModified(std::string filename);
bool fileExist(std::string filename);
int newStep();
void appendOption(std::string &str, const std::string &what, const int val);
void appendOption(std::string &str, const std::string &what);
void GmshDisplay(onelab::remoteNetworkClient *loader, std::string fileName, std::vector<std::string> choices);
void GmshDisplay(onelab::remoteNetworkClient *loader, std::string modelName, std::string fileName);
std::string getCurrentWorkdir();
std::string getUserHomedir();

int systemCall(std::string cmd);

typedef std::vector <std::vector <double> > array;
array read_array(std::string filename, char sep);
double find_in_array(const int i, const int j, const std::vector <std::vector <double> > &data);
std::vector<double> extract_column(const int j, array data);

class MetaModel : public onelab::localClient {
private:
  std::vector<onelab::client *> _clients;
  void registerClients();
 public:
 MetaModel(const std::string &commandLine, const std::string &cname, const std::string &fname, const int number) 
   : localClient(commandLine) {
    clientName = cname;
    genericNameFromArgs = fname;
    modelNumberFromArgs = number;
    registerClients();
  }
  ~MetaModel(){}

  std::string genericNameFromArgs, clientName;
  int modelNumberFromArgs;
  void registerClient(onelab::client *pName); 
  void initialize();
  void initializeClients();
  void analyze(); // the following 2 functions are defined by the user in a separate file
  void compute();
};

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
  //std::vector<std::string> getChoices(const std::string paramName);
  std::string stateToChar();  
  std::string showParamSpace();
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

class InterfacedClient : public onelab::localClient { 
  // n'utilise pas localNetworkClient::run
  // n'a donc pas pas besoin de _initializeCommand, _analyzeCommand, _computeCommand
  // les options sont transférées à compute() sans passer par le serveur
private:
  std::string _commandLine, _fileName, _extension, _options;
  std::set<std::string, ShortNameLessThan> _parameters;
  bool analyze_oneline(std::string line, std::ifstream &infile) ;
  bool analyze_ifstatement(std::ifstream &infile, bool condition) ;
  bool analyze_onefile(std::string ifilename);
  bool convert_oneline(std::string line, std::ifstream &infile, std::ofstream &outfile);
  bool convert_ifstatement(std::ifstream &infile, std::ofstream &outfile, bool condition) ;
  bool convert_onefile(std::string ifileName, std::ofstream &outfile);
public:
 InterfacedClient(const std::string &name, const std::string &commandLine, const std::string &extension) 
   : onelab::localClient(name), _commandLine(commandLine), _extension(extension) {}
  ~InterfacedClient(){}
  void setCommandLine(const std::string &cmd){ _commandLine.assign(cmd); }
  //std::string getFileName();
  void setFileName(const std::string &nam);
  void setLineOptions(const std::string &opt) { _options.assign(opt); }

  void initialize(); 
  void analyze();
  void convert();
  void compute();
};

class EncapsulatedClient : public onelab::localNetworkClient { 
  // utilise localNetworkClient::run
private:
  std::string _extension;
public:
 EncapsulatedClient(const std::string &name, const std::string &commandLine,  const std::string &extension) 
   : onelab::localNetworkClient(name,commandLine), _extension(extension) {}
  ~EncapsulatedClient(){}
  void setExtension(const std::string ext) { _extension.assign(ext); }
  std::string getExtension() { return _extension; }
  void setFileName(const std::string &nam); 
  std::string getFileName();
  void setLineOptions(const std::string &opt);

  void initialize(); 
  void analyze();
  void compute() ;
};

#endif



/* class InterfacedElmer : public InterfacedClient { */
/* public: */
/*   InterfacedElmer(const std::string &name) : InterfacedClient(name) { */
/*     setExtension(".sif"); */
/*     setCommandLine("ElmerSolver"); */
/*   } */
/*   ~InterfacedElmer(){} */
/* }; */

/* class LoadedMetaModel : public EncapsulatedClient { */
/* public: */
/*  LoadedMetaModel(const std::string &commandLine, const std::string fileName, const int modelNumber) */
/*    : EncapsulatedClient("metamodel",commandLine) { */
/*     setExtension(""); */
/*     setInitializeCommand("-h"); */
/*     setCheckCommand("-a"); */
/*     setComputeCommand(""); */
/*   } */
/*   ~LoadedMetaModel(){} */
/* }; */

/* class InterfacedElast : public InterfacedClient { */
/* public: */
/*   InterfacedElast(const std::string &name) : InterfacedClient(name) { */
/*     setExtension(".dat"); */
/*     setCommandLine("$GMSH_DIR/utils/api_demos/build/mainElasticity"); */
/*   } */
/*   ~InterfacedElast(){} */
/* }; */

/* class EncapsulatedGmsh : public EncapsulatedClient { */
/* public: */
/*   EncapsulatedGmsh(const std::string &name) : EncapsulatedClient(name,"gmsh") { */
/*     setExtension(".geo"); */
/*     // setInitializeCommand(""); */
/*     // setCheckCommand("-"); */
/*     // setComputeCommand("-3"); */
/*   } */
/*   ~EncapsulatedGmsh(){} */
/* }; */
