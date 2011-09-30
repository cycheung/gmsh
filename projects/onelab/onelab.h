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
    parameter(const std::string &name, const std::string &shortHelp="", 
              const std::string &help="")
      : _name(name), _shortHelp(shortHelp), _help(help){}
    void setShortHelp(std::string &shortHelp){ _shortHelp = shortHelp; }
    void setHelp(std::string &help){ _help = help; }
    void setClients(std::set<std::string> &clients){ _clients = clients; }
    void addClient(std::string &client){ _clients.insert(client); }
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
    char charSep(){ return '\0'; }
    std::string sanitize(const std::string &in)
    {
      std::string out(in);
      for(unsigned int i = 0; i < in.size(); i++)
        if(out[i] == charSep()) out[i] = ' ';
      return out;
    }
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
      : parameter(name), _value(0.), _defaultValue(0.), 
        _min(1.), _max(0.), _step(0.)
    {
    }
    number(const std::string &name, double defaultValue,
           const std::string &shortHelp="", const std::string &help="") 
      : parameter(name, shortHelp, help), _value(defaultValue), 
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
      sstream << getType() << charSep() << sanitize(getName()) << charSep() 
              << sanitize(getShortHelp()) << charSep() 
              << sanitize(getHelp()) << charSep() << _value << charSep()
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
      : parameter(name), _value(""), _defaultValue("")
    {
    }
    string(const std::string &name, const std::string &defaultValue, 
           const std::string &shortHelp="", const std::string &help="") 
      : parameter(name, shortHelp, help), _value(defaultValue), 
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
      sstream << getType() << charSep() << sanitize(getName()) << charSep() 
              << sanitize(getShortHelp()) << charSep() 
              << sanitize(getHelp()) << charSep() << sanitize(_value) << charSep()
              << sanitize(_defaultValue) << charSep()
              << _choices.size() << charSep();
      for(unsigned int i = 0; i < _choices.size(); i++)
        sstream << sanitize(_choices[i]) << charSep();
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
      : parameter(name), _value(""), _defaultValue("")
    {
    }
    region(const std::string &name, const std::string &defaultValue,
           const std::string &shortHelp="", const std::string &help="") 
      : parameter(name, shortHelp, help), _value(defaultValue), 
        _defaultValue(defaultValue)
    {
    }
    parameterType getType() const { return REGION; }
    const std::string &getValue() const { return _value; }
    const std::string &getDefaultValue() const { return _defaultValue; }
    std::string toChar()
    {
      std::ostringstream sstream;
      sstream << getType() << charSep() << sanitize(getName()) << charSep() 
              << sanitize(getShortHelp()) << charSep() << sanitize(getHelp()) 
              << charSep() << _value << charSep() << _defaultValue << charSep()
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
      : parameter(name), _value(""), _defaultValue("")
    {
    }
    function(const std::string &name, const std::string &defaultValue,
             const std::string &shortHelp="", const std::string &help="") 
      : parameter(name, shortHelp, help), _value(defaultValue), 
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
      sstream << getType() << charSep() << sanitize(getName()) << charSep() 
              << sanitize(getShortHelp()) << charSep() 
              << sanitize(getHelp()) << charSep() << sanitize(_value) << charSep()
              << sanitize(_defaultValue) << charSep()
              << _pieceWiseValues.size() << charSep();
        for(std::map<std::string, std::string>::const_iterator it =
              _pieceWiseValues.begin(); it != _pieceWiseValues.end(); it++)
          sstream << sanitize(it->first) << charSep() 
                  << sanitize(it->second) << charSep();
      sstream << _choices.size() << charSep();
      for(unsigned int i = 0; i < _choices.size(); i++)
        sstream << sanitize(_choices[i]) << charSep();
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
    template <class T> bool _set(T &p, std::set<T*, parameterLessThan> &parameters)
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
      return true;
    }
    // get the parameter matching the given name, or all the
    // parameters in the category if no name is given
    template <class T> bool _get(std::vector<T> &p, const std::string &name,
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
      return true;
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
    bool set(number &p){ return _set(p, _numbers); }
    bool set(string &p){ return _set(p, _strings); }
    bool set(region &p){ return _set(p, _regions); }
    bool set(function &p){ return _set(p, _functions); }
    bool get(std::vector<number> &p, const std::string &name="")
    {
      return _get(p, name, _numbers); 
    }
    bool get(std::vector<string> &p, const std::string &name="")
    {
      return _get(p, name, _strings); 
    }
    bool get(std::vector<region> &p, const std::string &name="")
    {
      return _get(p, name, _regions); 
    }
    bool get(std::vector<function> &p, const std::string &name="")
    {
      return _get(p, name, _functions);
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
    // the connected clients, indexed by name
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
    template <class T> bool set(T &p)
    {
      // this needs to be locked to avoid race conditions when several
      // clients try to set a parameter at the same time
      return _parameterSpace.set(p); 
    }
    template <class T> bool get(std::vector<T> &p, const std::string &name="")
    {
      return _parameterSpace.get(p, name); 
    }
    bool registerClient(const std::string &name, client *c)
    {
      if(_clients.count(name)) return false;
      _clients[name] = c;
      return true;
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
  public:
    client(const std::string &name) : _name(name){}
    virtual ~client(){}
    virtual bool analyze(){ return false; }
    virtual bool compute(const std::string &what){ return false; }
  };

  class localClient : public client{
  protected:
    // the pointer to the server
    server *_server;
  public:
    localClient(const std::string &name) : client(name)
    {
      _server = server::instance();
      _server->registerClient(_name, this);
    }
    virtual ~localClient(){}
    template <class T> bool set(T &parameter)
    {
      parameter.addClient(_name);
      _server->set(parameter);
      return true;
    }
    template <class T> bool get(std::vector<T> &parameters,
                                const std::string &name="")
    {
      _server->get(parameters, name);
      return true;
    }
  };

  class localNetworkClient : public localClient{
  private:
    std::string _commandLine;
  public:
    localNetworkClient(const std::string &name, const std::string &commandLine)
      : localClient(name), _commandLine(commandLine)
    {
      // new GmshServer()
    }
    virtual ~localNetworkClient()
    {
      // delete GmshServer
    }
    virtual bool analyze()
    {
      // launch distant client using command line and the server address
      // send analyze request
      // while(client stop)
      //   wait for data from client
    }
    virtual bool compute(const std::string &what)
    {
      // launch remote client using command line and the server address
      // send compute request
      // while(client stop)
      //   wait for data from client
    }
  };

  class remoteNetworkClient : public client{
  private:
    std::string _serverAddress;
  public:
    remoteNetworkClient(const std::string &name, const std::string &serverAddress)
      : client(name), _serverAddress(serverAddress)
    {
        /*
        _gmshClient = new GmshClient();
        if(_gmshClient->Connect(_serverAddress.c_str()) < 0){
          throw "Could not connect";
          delete _gmshClient;
          _gmshClient = 0;
        }
        else{
          _gmshClient->Start();
        }
        */
    }
    virtual ~remoteNetworkClient()
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
    template <class T> bool set(T &parameter)
    {
      parameter.addClient(_name);
      // send over socket
      return false;
    }
    template <class T> bool get(std::vector<T> &parameters, 
                                const std::string &name="")
    {
      // ask for data
      // while(1){
      //  wait for answer
      // }
      return false;
    }
  };

}
