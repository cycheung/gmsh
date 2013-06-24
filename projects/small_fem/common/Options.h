#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <string>
#include <vector>
#include <map>

/**
   @class Options
   @brief Handel options

   An Option is a pair composed two strings:
   @li The first one is called the @em option
   @li The second one is called the @em value

   This class can store all the pairs (option, value) for latter access.
*/

class Options{
 private:
  std::multimap<std::string, std::string>* optionMap;

 public:
  Options(size_t nArg, const char *const *const arg);
  ~Options(void);

  std::vector<std::string> getValue(std::string option) const;
  std::string toString(void) const;

  static void cStyle(std::vector<std::string>& vec,
                     char** vecCStyle,
                     size_t offset);
};

/**
   @fn Options::Options
   @param arg A vector of string
   @param nArg The size of this vector

   Instanciates a new Option with the given vector of string.

   The options will taken in the folowing way:
   (option, value)[i] = (arg[2 * i] + arg[2 * i + 1]) for i = {0, ..., nArg / 2}.

   All the paires (option, value) will be stored for latter use.

   @see Options::getValue
   **

   @fn Options::~Options
   Deletes this Options
   **

   @fn Options::getValue
   @param string A string used as an option
   @return
   Returns a vector with the values stored for the given option.@n
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

   @warning
   The char* are @em bounded to the given vector.
   **
 */

#endif
