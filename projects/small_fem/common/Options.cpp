#include <sstream>
#include "Exception.h"
#include "Options.h"

using namespace std;

Options::Options(size_t nArg, const char *const *const arg){
  // Get Size
  const size_t size = nArg / 2;

  // Init Map
  optionMap = new multimap<string, string>;

  // Fill map
  for(size_t i = 0; i < size; i++)
    optionMap->insert(pair<string, string>(string(arg[2 * i]),
                                           string(arg[2 * i + 1])));
}

Options::~Options(void){
  delete optionMap;
}

vector<string> Options::getValue(string option) const{
  // Get iterators to option in optionMap
  const pair<multimap<string, string>::const_iterator,
             multimap<string, string>::const_iterator> itRef =
    optionMap->equal_range(option);

  multimap<string, string>::const_iterator it;

  // Get number of element in range
  size_t N = 0;

  for(it = itRef.first; it != itRef.second; it++)
    N++;

  // Reset 'it'
  it = itRef.first;

  // Init Vector to Return
  vector<string> v(N);

  for(size_t i = 0; it != itRef.second; it++, i++)
    v[i] = it->second;

  // Return
  return v;
}

void Options::cStyle(vector<string>& vec,
                     char** vecCStyle,
                     size_t offset){
  // Get size
  const size_t N = vec.size();

  for(size_t i = 0; i < N; i++)
    vecCStyle[i + offset] = const_cast<char*>(vec[i].c_str());
}

string Options::toString(void) const{
  stringstream stream;

  const multimap<string, string>::const_iterator end = optionMap->end();
  multimap<string, string>::const_iterator        it = optionMap->begin();

  for(; it != end; it++)
    stream << "("
           << it->first
           << ", "
           << it->second
           << ")"
           << endl;

  return stream.str();
}
