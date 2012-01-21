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

int getOptions(int argc, char *argv[],int &modelNumber, bool &analyzeOnly, std::string &name, std::string &sockName);
std::string itoa(const int i);
int metamodel(int modelNumber);
int analyze();
int compute();
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
  std::vector<std::string> getChoices(const std::string paramName);
  std::string stateToChar();  
  std::string showParamSpace();
  bool menu(std::string options, std::string modelName);
};


class InterfacedClient : public onelab::localClient {
private:
  std::string _name, _extension, _commandLine;
public:
  InterfacedClient(const std::string &name) : onelab::localClient(name) {}
  ~InterfacedClient(){}
  void setName(const std::string nam) { _name.assign(nam); }
  void setExtension(const std::string ext) { _extension.assign(ext); }
  void setCommandLine(const std::string cmd) { _commandLine.assign(cmd); }

  bool analyze(const std::string options, const std::string modelName);
  bool analyze_oneline(std::string line, std::ifstream &infile) ;
  bool analyze_ifstatement(std::ifstream &infile, bool condition) ;
  bool analyze_onefile(std::string ifilename);
  bool convert(std::string modelName);
  bool convert_oneline(std::string line, std::ifstream &infile, std::ofstream &outfile);
  bool convert_ifstatement(std::ifstream &infile, std::ofstream &outfile, bool condition) ;
  bool convert_onefile(std::string ifileName, std::ofstream &outfile);

  bool run(const std::string options, const std::string modelName) ;
};

class InterfacedElmer : public InterfacedClient {
public:
  InterfacedElmer(const std::string &name) : InterfacedClient(name) {
    setName(name);
    setExtension(".sif");
    setCommandLine("ElmerSolver");
  }
  ~InterfacedElmer(){}
};

class InterfacedGmsh : public InterfacedClient {
 public:
 InterfacedGmsh(const std::string &name) : InterfacedClient(name) {
    setName(name);
    setExtension(".geo");
    setCommandLine("gmsh");
  }
  ~InterfacedGmsh(){}
};

class InterfacedGetdp : public InterfacedClient {
public:
  InterfacedGetdp(const std::string &name) : InterfacedClient(name) {
    setName(name);
    setExtension(".pro");
    setCommandLine("getdp");
  }
  ~InterfacedGetdp(){}
};

class InterfacedElast : public InterfacedClient {
public:
  InterfacedElast(const std::string &name) : InterfacedClient(name) {
    setName(name);
    setExtension("");
    setCommandLine("$GMSH_DIR/utils/api_demos/build/mainElasticity");
  }
  ~InterfacedElast(){}
};

class EncapsulatedTest : public onelab::localNetworkClient {
public:
  EncapsulatedTest(const std::string &name) : onelab::localNetworkClient(name,name) {}
  ~EncapsulatedTest(){}
  bool analyze(const std::string options, const std::string modelName);
  bool run(const std::string options, const std::string modelName);
};


class EncapsulatedGmsh : public onelab::localNetworkClient {
public:
  EncapsulatedGmsh(const std::string &name) : onelab::localNetworkClient(name,"gmsh") {}
  ~EncapsulatedGmsh(){}
  int  getVerbosity();
  bool analyze(const std::string options, const std::string modelName);
  bool run(const std::string options, const std::string modelName);
};


class EncapsulatedGetdp : public onelab::localNetworkClient {
public:
  EncapsulatedGetdp(const std::string &name) : onelab::localNetworkClient(name,"getdp") {}
  ~EncapsulatedGetdp(){}
  int  getVerbosity();
  bool analyze(const std::string options, const std::string modelName);
  bool run(const std::string options, const std::string modelName);
};

class MetaModel : public onelab::localNetworkClient {
 public:
  MetaModel(const std::string &name) : onelab::localNetworkClient(name,name) {}
  ~MetaModel(){}
  bool analyze(const std::string options, const std::string modelName);
  bool run(const std::string options, const std::string modelName);
};

#endif
