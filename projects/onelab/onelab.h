// ONELAB - Copyright (C) 2011 ULg/UCL
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, and/or sell copies of the
// Software, and to permit persons to whom the Software is furnished
// to do so, provided that the above copyright notice(s) and this
// permission notice appear in all copies of the Software and that
// both the above copyright notice(s) and this permission notice
// appear in supporting documentation.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR
// ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY
// DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
// WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
// ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
// OF THIS SOFTWARE.

#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <sstream>
//#include "GmshSocket.h"

namespace onelab{

  typedef enum { NUMBER = 1, STRING = 2, REGION = 3, FUNCTION = 4 } parameterType;

  class number;
  class string;
  class region;
  class function;
  class client;
  class server;

  // The base parameter class.
  class parameter{
  private:
    // the name of the parameter, including its "path" in the
    // parameter hierarchy. The path separator '/' can be followed by
    // a number to force ordering (hence a parameter name cannot start
    // with a number).
    std::string _name;
    // optional help strings
    std::string _shortHelp, _help;
    // client code(s) for which this parameter makes sense
    std::set<std::string> _clients;
  public:
    parameter(const std::string &client, const std::string &name, 
              const std::string &shortHelp="", const std::string &help="")
      : _name(name), _shortHelp(shortHelp), _help(help)
    {
      _clients.insert(client);
    }
    void setShortHelp(std::string &shortHelp){ _shortHelp = shortHelp; }
    void setHelp(std::string &help){ _help = help; }
    void setClients(std::set<std::string> &clients){ _clients = clients; }
    void addClients(std::set<std::string> &clients)
    { 
      _clients.insert(clients.begin(), clients.end()); 
    }
    virtual parameterType getType() const = 0;
    std::string getTypeAsString()
    {
      std::ostringstream sstream;
      sstream << getType();
      return sstream.str();
    }
    const std::string &getName() const { return _name; }
    const std::string &getShortHelp() const { return _shortHelp; }
    const std::string &getHelp() const { return _help; }
    const std::set<std::string> &getClients() const { return _clients; }
    std::string charSep(){ return "|"; }
    virtual std::string toChar() = 0;
    virtual void fromChar(const std::string &c){}
  };
  
  class parameterLessThan{
  public:
    bool operator()(const parameter *p1, const parameter *p2) const
    {
      return p1->getName() < p2->getName();
    }
  };

  // The number class. Numbers are stored internally as double
  // precision real numbers. Currently all more complicated types
  // (complex numbers, vectors, etc.) are supposed to be encapsulated
  // in functions. We will probably add more base types in the future
  // to make the interface nicer.
  class number : public parameter{
  private:
    double _value;
    double _defaultValue, _min, _max, _step;
    std::vector<double> _choices;
  public:
    number(const std::string &name) 
      : parameter("", name), _value(0.), _defaultValue(0.), 
        _min(1.), _max(0.), _step(0.)
    {
    }
    number(const std::string &client, const std::string &name, double defaultValue,
           const std::string &shortHelp="", const std::string &help="") 
      : parameter(client, name, shortHelp, help), _value(defaultValue), 
        _defaultValue(defaultValue), _min(1.), _max(0.), _step(0.)
    {
    }
    void setValue(double value){ _value = value; }
    void setMin(double min){ _min = min; }
    void setMax(double max){ _min = max; }
    void setStep(double step){ _step = step; }
    void setChoices(std::vector<double> &choices){ _choices = choices; }
    parameterType getType() const { return NUMBER; }
    double getValue() const { return _value; }
    double getDefaultValue() const { return _defaultValue; }
    std::string toChar()
    {
      std::ostringstream sstream;
      sstream << getType() << charSep() << getName() << charSep() 
              << getShortHelp() << charSep() << getHelp() << charSep()
              << _value << charSep()
              << _defaultValue << charSep()
              << _min << charSep() << _max << charSep() << _step << charSep()
              << _choices.size() << charSep();
      for(unsigned int i = 0; i < _choices.size(); i++)
        sstream << _choices[i] << charSep();
      return sstream.str();
    }
  };

  // The string class.
  class string : public parameter{
  private:
    std::string _value;
    std::string _defaultValue;
    std::vector<std::string> _choices;
  public:
    string(const std::string &name) 
      : parameter("", name), _value(""), _defaultValue("")
    {
    }
    string(const std::string &client, const std::string &name, 
           const std::string &defaultValue, const std::string &shortHelp="",
           const std::string &help="") 
      : parameter(client, name, shortHelp, help), _value(defaultValue), 
        _defaultValue(defaultValue)
    {
    }
    void setValue(const std::string &value){ _value = value; }
    void setChoices(std::vector<std::string> &choices){ _choices = choices; }
    parameterType getType() const { return STRING; }
    const std::string &getValue() const { return _value; }
    const std::string &getDefaultValue() const { return _defaultValue; }
    std::string toChar()
    {
      std::ostringstream sstream;
      sstream << getType() << charSep() << getName() << charSep() 
              << getShortHelp() << charSep() << getHelp() << charSep()
              << _value << charSep()
              << _defaultValue << charSep()
              << _choices.size() << charSep();
      for(unsigned int i = 0; i < _choices.size(); i++)
        sstream << _choices[i] << charSep();
      return sstream.str();
    }
  };

  // The region class. A region can be any kind of geometrical entity.
  class region : public parameter{
  private:
    std::string _value, _defaultValue;
    std::vector<std::string> _choices;
  public:
    region(const std::string &name) 
      : parameter("", name), _value(""), _defaultValue("")
    {
    }
    region(const std::string &client, const std::string &name, 
           const std::string &defaultValue, const std::string &shortHelp="",
           const std::string &help="") 
      : parameter(client, name, shortHelp, help), _value(defaultValue), 
        _defaultValue(defaultValue)
    {
    }
    parameterType getType() const { return REGION; }
    const std::string &getValue() const { return _value; }
    const std::string &getDefaultValue() const { return _defaultValue; }
    std::string toChar()
    {
      std::ostringstream sstream;
      sstream << getType() << charSep() << getName() << charSep() 
              << getShortHelp() << charSep() << getHelp() << charSep()
              << _value << charSep()
              << _defaultValue << charSep()
              << _choices.size() << charSep();
      for(unsigned int i = 0; i < _choices.size(); i++)
        sstream << _choices[i] << charSep();
      return sstream.str();
    }
  };

  // The (possibly piece-wise defined on regions) function
  // class. Currently functions are entirely client-dependent: they
  // are just represented internally as strings. Again, we might want
  // to specialize in the future to make the interface more refined.
  class function : public parameter{
  private:
    std::string _value, _defaultValue;
    std::map<std::string, std::string> _pieceWiseValues;
    std::vector<std::string> _choices;
  public:
    function(const std::string &name) 
      : parameter("", name), _value(""), _defaultValue("")
    {
    }
    function(const std::string &client, const std::string &name, 
             const std::string &defaultValue, const std::string &shortHelp="",
             const std::string &help="") 
      : parameter(client, name, shortHelp, help), _value(defaultValue), 
        _defaultValue(defaultValue)
    {
    }
    void setValue(const std::string &value, const std::string &region="")
    {
      if(region.empty())
        _value = value;
      else
        _pieceWiseValues[region] = value;
    }
    parameterType getType() const { return FUNCTION; }
    const std::string getValue(const std::string &region="") const
    {
      if(region.size()){
        std::map<std::string, std::string>::const_iterator it =
          _pieceWiseValues.find(region);
        if(it != _pieceWiseValues.end()) return it->second;
        return "";
      }
      else return _value; 
    }
    const std::map<std::string, std::string> &getPieceWiseValues() const 
    {
      return _pieceWiseValues; 
    }
    const std::string &getDefaultValue() const { return _defaultValue; }
    std::string toChar()
    {
      std::ostringstream sstream;
      sstream << getType() << charSep() << getName() << charSep() 
              << getShortHelp() << charSep() << getHelp() << charSep()
              << _value << charSep()
              << _defaultValue << charSep()
              << _pieceWiseValues.size() << charSep();
        for(std::map<std::string, std::string>::const_iterator it =
              _pieceWiseValues.begin(); it != _pieceWiseValues.end(); it++)
          sstream << it->first << charSep() << it->second << charSep();
      sstream << _choices.size() << charSep();
      for(unsigned int i = 0; i < _choices.size(); i++)
        sstream << _choices[i] << charSep();
      return sstream.str();
    }
  };

  // The parameter space, i.e., the set of parameters stored and
  // handled by the onelab server.
  class parameterSpace{
  private:
    std::set<number*, parameterLessThan> _numbers;
    std::set<string*, parameterLessThan> _strings;
    std::set<region*, parameterLessThan> _regions;
    std::set<function*, parameterLessThan> _functions;
    // set a parameter in the parameter space; if it already exists,
    // use the new value but make sure to add new clients if necessary
    template <class T> void _set(T &p, std::set<T*, parameterLessThan> &parameters)
    {
      std::set<std::string> clients;
      typename std::set<T*, parameterLessThan>::iterator it = parameters.find(&p);
      if(it != parameters.end()){
        parameters.erase(it);
        T* oldp = *it;
        clients = oldp->getClients();
        delete oldp;
      }
      T* newp = new T(p);
      newp->addClients(clients);
      parameters.insert(newp);
    }
    // get the parameter matching the given name, or all the
    // parameters in the category if no name is given
    template <class T> void _get(std::vector<T> &p, const std::string &name,
                                 std::set<T*, parameterLessThan> &parameters)
    {
      if(name.empty()){
        for(typename std::set<T*, parameterLessThan>::iterator it = parameters.begin();
            it != parameters.end(); it++)
          p.push_back(**it);
      }
      else{
        T tmp(name);
        typename std::set<T*, parameterLessThan>::iterator it = parameters.find(&tmp);
        if(it != parameters.end())
          p.push_back(**it);
      }
    }
  public:
    parameterSpace(){}
    ~parameterSpace()
    {
      for(std::set<number*, parameterLessThan>::iterator it = _numbers.begin();
          it != _numbers.end(); it++)
        delete *it;
      for(std::set<string*, parameterLessThan>::iterator it = _strings.begin();
          it != _strings.end(); it++)
        delete *it;
      for(std::set<region*, parameterLessThan>::iterator it = _regions.begin();
          it != _regions.end(); it++)
        delete *it;
      for(std::set<function*, parameterLessThan>::iterator it = _functions.begin();
          it != _functions.end(); it++)
        delete *it;
    }
    void set(number &p){ _set(p, _numbers); }
    void set(string &p){ _set(p, _strings); }
    void set(region &p){ _set(p, _regions); }
    void set(function &p){ _set(p, _functions); }
    void get(std::vector<number> &p, const std::string &name="")
    {
      _get(p, name, _numbers); 
    }
    void get(std::vector<string> &p, const std::string &name="")
    {
      _get(p, name, _strings); 
    }
    void get(std::vector<region> &p, const std::string &name="")
    {
      _get(p, name, _regions); 
    }
    void get(std::vector<function> &p, const std::string &name="")
    {
      _get(p, name, _functions);
    }
  };

  // The onelab server: a singleton that stores the parameter space
  // and interacts with onelab clients.
  class server{
  private:
    // the unique server
    static server *_server;
    // the address of the server
    std::string _address;
    // the connected clients, indexed by address
    std::map<std::string, client*> _clients;
    // the parameter space
    parameterSpace _parameterSpace;
  public:
    server(const std::string &address="") : _address(address) {}
    ~server(){}
    static server *instance(const std::string &address="")
    {
      if(!_server) _server = new server(address);
      return _server;
    }
    template <class T> void set(T &p){ _parameterSpace.set(p); }
    template <class T> void get(std::vector<T> &p, const std::string &name="")
    {
      _parameterSpace.get(p, name); 
    }
    void registerClient(const std::string &name, client *c)
    {
      _clients[name] = c;
    }
    typedef std::map<std::string, client*>::iterator citer;
    citer firstClient(){ return _clients.begin(); }
    citer lastClient(){ return _clients.end(); }
    citer findClient(const std::string &name){ return _clients.find(name); }
  };
    
  // The onelab client: a class that communicates with the onelab
  // server. Each client should be derived from this one.
  class client{
  protected:
    // the name of the client
    std::string _name;
    // the command-line to run a separate code if necessary
    std::string _commandLine;
    // the GmshClient if the client is distant
    //GmshClient *_gmshClient;
    // the pointer to the server if the client is local
    server *_server;
  public:
    client(const std::string &name, const std::string &commandLine="",
           const std::string &serverAddress="") 
      : _name(name), _commandLine(commandLine), /*_gmshClient(0),*/ _server(0)
    {
      if(serverAddress.empty()){
        // client and server are in the same process
        _server = server::instance();
        _server->registerClient(_name, this);
      }
      else{
        /*
        _gmshClient = new GmshClient();
        if(_gmshClient->Connect(serverAddress.c_str()) < 0){
          throw "Could not connect";
          delete _gmshClient;
          _gmshClient = 0;
        }
        else{
          _gmshClient->Start();
        }
        */
      }
    }
    virtual ~client()
    {
      /*
      if(_gmshClient){
        _gmshClient->Stop();
        _gmshClient->Disconnect();
        delete _gmshClient;
        _gmshClient = 0;
      }
      */
    }
    // run client code
    virtual void run(){}
    template <class T> bool send(T &parameter)
    {
      if(_server){
        _server->set(parameter);
        return true;
      }
      else{
        // send over socket
        return false;
      }
    }
    template <class T> bool receive(std::vector<T> &parameters, 
                                    const std::string &name="")
    {
      if(_server){
        _server->get(parameters, name);
        return true;
      }
      else{
        // receive over socket
        return false;
      }
    }
  };

  // example of a command-line client, that communicates with a
  // remote code using GmshSocket
  class commandLineClient : public client {
  public:
    commandLineClient(const std::string &name, const std::string &commandLine) 
      : client(name, commandLine, "")
    {
    }
    virtual void run()
    {
      // new GmshServer()
      // GmshServer->Start(_commandLine)
      // do{
      //   listen to client and respond to its queries
      // } while(we get stop signal from client)
      // GmshServer->Shutdown()
    }
  };

  // example of a local client, that can interact directly with the
  // parameters stored in the server
  class gmshClient : public client {
  public:
    gmshClient() : client("gmsh", "", "") {}
    virtual void run()
    {
      // get parameters from _server
      // export .geo.opt
      // e.g. reparse, mesh, saveMesh for current model
    }
  };

  // example of a distant client, that interacts with the
  // commandLineClient through GmshSocket
  class getdpClient : public client {
  public:
    getdpClient(const std::string &serverAddress) 
      : client("getdp", "", serverAddress) {}
  };

}
