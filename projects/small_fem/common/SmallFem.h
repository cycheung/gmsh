#ifndef _SMALLFEM_H_
#define _SMALLFEM_H_

#include "Options.h"

/**
   @class SmallFem
   @brief SmallFem Initialize and finalize

   Initialize and Finalize SmallFem.
   It also gives access to Options.
*/

class SmallFem{
 private:
  static bool     initOne;
  static bool     finaOne;
  static Options* option;

 public:
   SmallFem(void);
  ~SmallFem(void);

  static void Initialize(int argc, char** argv);
  static void Finalize(void);

  static Options& getOptions(void);
};

/**
   @fn SmallFem::SmallFem
   Instantiates a new SmallFem

   Not needed, since SmallFem is bunch of static class
   **

   @fn SmallFem::~SmallFem
   Deletes this SmallFem

   Not needed, since SmallFem is bunch of static class
   **

   @fn SmallFem::Initialize
   @param argv A vector of char*
   @param argc The size of the previous vector

   Class method initializes SmallFem. Moreover it
   takes {argv[1], ..., argv[argc - 1]} as Options.

   @see SmallFem::getOption()
   **

   @fn SmallFem::Finalize
   Class method finalizes SmallFem
   **

   @fn SmallFem::getOptions
   @return Returns the Options given in SmallFem::Initialize()
   **
*/

#endif
