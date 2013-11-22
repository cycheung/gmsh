#include <sstream>
#include <set>
#include <list>

#include "Options.h"

using namespace std;

Options::Options(int argc,char** argv, const string& keywords){
  // Set of keywords //
  size_t      it  = 0;
  size_t      pos = 0;
  set<string> key;

  do{
    it = keywords.find(",", pos);

    if(it != string::npos)
      key.insert(keywords.substr(pos, it - pos));

    else
      key.insert(keywords.substr(pos, it));

    pos = it + 1;
  }while(it != string::npos);

  // Go to first option //
  int offset = 0;
  while(offset < argc && key.count(string(argv[offset])) != 1)
    offset++;

  // Parse options if any //
  list<pair<string, string> > list;
  string opt;

  if(offset != argc){
    opt = string(argv[offset]);

    for(int i = offset + 1; i < argc; i++){
      if(key.count(string(argv[i])) == 1)
        // A new option is found
        opt = string(argv[i]);

      else
        // Enlist this option arguments
        list.push_back(pair<string, string>(opt, string(argv[i])));
    }
  }

  // Save fifo
  optionMap = new multimap<string, string>(list.begin(), list.end());
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
