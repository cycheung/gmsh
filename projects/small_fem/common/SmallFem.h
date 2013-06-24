#ifndef _SMALLFEM_H_
#define _SMALLFEM_H_

#include "Options.h"

/**
   @class SmallFem
   @brief SmallFem Initialize and finalize

   SmallFem Initialize and Finalize
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
   Class method initializing SmallFem
   **

   @fn SmallFem::Finalize
   Class method finalizing SmallFem
*/

#endif
