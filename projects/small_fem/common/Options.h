#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <string>
#include <vector>
#include <map>

/**
   @class Options
   @brief Handel options

   An Option is a pair composed two strings:
   @li The first one is called the option
   @li The second one is called the value

   This class can store all the pairs (option, value) for latter access.

   The possible values of the options are given by a set of keywords.
   These keywords are represented by a comma separated string.
*/

class Options{
 private:
  std::multimap<std::string, std::string>* optionMap;

 public:
  Options(int argc, char** argv, const std::string& keywords);
  ~Options(void);

  std::vector<std::string> getValue(std::string option) const;
  std::string toString(void) const;

  static void cStyle(std::vector<std::string>& vec,
                     char** vecCStyle,
                     size_t offset);
};

/**
   @fn Options::Options
   @param argv A vector of char*
   @param argc The size of the previous vector
   @param keywords A string

   Instanciates a new Options.

   This method will parse the given char* and will store all the
   (option, value) pairs, with the given keywords as possible options.

   All the paires (option, value) will be stored for latter use.

   @see Options::getValue
   **

   @fn Options::~Options
   Deletes this Options
   **

   @fn Options::getValue
   @param option A string used as an option
   @return
   Returns a vector with the values stored for the given option.
   Those values have been stored during instanciation of this Options.

   @see Options::Options
   **

   @fn Options::toString
   @return
   Returns a string with the stored options
   **

   @fn Options::cStyle
   @param vec A vector of string
   @param vecCStyle An allocated array of char* of size at least vec.size()
   @param offset An positive integer

   Convert the given vector of string into a char**.
   The given char** will populated with the char* corresponding
   to each entry of the given vector.
   vecCStyle will be populated starting at offset.

   The char* are bounded to the given vector.
   **
 */

#endif
