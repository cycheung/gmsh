#include "onelab.h"

onelab::server *onelab::server::_server = 0;

int main()
{
  // server starts client, which sends its parameters
  onelab::number mu0("getdp", "/Physical parameters/mu0", 2.23,
                       "Magnetic permeability of vaccum");
  onelab::number sigma("getdp", "/Geometry/numEncoches", 2.);
  onelab::number sigma2("gmsh", "/Geometry/numEncoches", 10.);
  onelab::region omega("getdp", "/Geometry/Omega", "100, 200, \"test\"");
  onelab::region gamma("getdp", "/Geometry/Gamma", "1000");
  onelab::function mu("getdp", "/Physical parameters/mu", "100 * mu0",
                      "Magnetic permeability");
  mu.setValue("200 * mu0", omega.getName());
  mu.setValue("300 * mu0", gamma.getName());
  onelab::client c("");//("localhost:4455");
  c.send(mu0);
  c.send(sigma);
  c.send(sigma2);
  c.send(omega);
  c.send(mu);
  
  // server changes parameters (e.g. graphically)

  // server starts client, client gets parameters and runs
  std::vector<onelab::number> numbers;
  onelab::client c2("");//("localhost:4455");
  c2.receive(numbers);
  for(unsigned int i = 0; i < numbers.size(); i++)
    std::cout << numbers[i].toChar() << std::endl;

}

#if 0

std::string next_token(const std::string &line, const std::string &delim, 
                       std::string::size_type &first)
{
  string::size_type last = line.find_first_of(delim, first);
  string next(line.substr(first, last));
  first = last;
  return next;
}

// server
main()
{
  onelab::server *s = new onelab::server();
  onelab::commandLineClient *c = new onelab::commandLineClient
    ("getdp", "getdp -onelab " + s->getAddress());
  serv->registerClient(c);

}

// client
main()
{
  // if argv == -onelab localhost:44266
  //   Msg::InitOnelab(localhost:44266)

  // Msg::InitOnelab(addr){
  //   _onelabClient = new remoteClient(addr);
  // }

  // in main:
  // _onelabClient->getParameters(&param)
  // if(param == parse)
  //   GetPDParse
  // _onelabClient->Stop()
  
  // in GetDPParse:
  //   Msg::SendParameter()


  onelab::getdpClient *c = new onelab::getdpClient(addr);
  c->run();

}

#endif
