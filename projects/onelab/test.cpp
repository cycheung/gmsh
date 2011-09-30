#include "onelab.h"

onelab::server *onelab::server::_server = 0;

int main()
{
  // server starts client, which sends its parameters
  onelab::number mu0("/Physical parameters/mu0", 2.23,
                       "Magnetic permeability of vaccum");
  onelab::number sigma("/Geometry/numEncoches", 2.);
  onelab::number sigma2("/Geometry/numEncoches", 10.);
  onelab::region omega("/Geometry/Omega", "100, 200, \"test\"");
  onelab::region gamma("/Geometry/Gamma", "1000");
  onelab::function mu("/Physical parameters/mu", "100 * mu0", "Magnetic permeability");
  mu.setValue("200 * mu0", omega.getName());
  mu.setValue("300 * mu0", gamma.getName());

  onelab::localClient c("getdp");
  c.set(mu0);
  c.set(sigma);
  c.set(sigma2);
  c.set(omega);
  c.set(mu);
  
  // server changes parameters (e.g. graphically)

  // server starts client, client gets parameters and runs
  onelab::localClient c2("getdp");
  std::vector<onelab::number> numbers;
  c2.get(numbers);
  for(unsigned int i = 0; i < numbers.size(); i++)
    std::cout << numbers[i].toChar() << std::endl;

  std::vector<onelab::function> functions;
  c2.get(functions);
  for(unsigned int i = 0; i < functions.size(); i++)
    std::cout << functions[i].toChar() << std::endl;

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
  onelab::localNetworkClient *c = new onelab::localNetworkClient
    ("getdp", "getdp -onelab " + s->getAddress());
  serv->registerClient(c);
}

// client
main()
{
  // if argv == -onelab localhost:44266
  // _onelabClient = new remoteNetworkClient("getdp", addr);

  // in parser, when we have a DefineVariable
  // if pre || cal || pos
  //   _onelabClient->get(parameters, name)
  // else
  //   _onelabClient->set(param)

  // delete _onelabClient;
}

#endif
