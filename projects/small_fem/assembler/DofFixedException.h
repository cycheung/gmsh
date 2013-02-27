#ifndef DOFFIXEDEXCEPTION_H_
#define DOFFIXEDEXCEPTION_H_

#include <exception>
#include "Dof.h"

/**
   @class DofFixedException
   @brief Exceptions indicating that a Dof has a fixed value

   Exceptions indicating that a Dof has a fixed value.@n
   The Dof and its value can be recovered by this exception.

   @warning
   DofFixedException%s are @em not Exception%s
 */

class DofFixedException: public std::exception{
 private:
  const Dof* dof;
  double     value;
  char*      string;

 public:
  DofFixedException(const Dof& dof, double value);

  virtual ~DofFixedException(void) throw();
  virtual const char* what(void) const throw();

  const Dof& getDof(void)   const throw();
  double     getValue(void) const throw();
};

/**
   @fn DofFixedException::DofFixedException
   @param dof A Dof
   @param value A real number

   Instantiates a new DofFixedException such that:
   @li DofFixedException::getDof() is the given Dof
   @li DofFixedException::getValue() is the given value
   **

   @fn DofFixedException::~DofFixedException
   Deletes this DofFixedException
   **

   @fn DofFixedException::what
   @return Returns the @em cause of the exception
   **

   @fn DofFixedException::getDof
   @return Returns the Dof given in the constructor
   **

   @fn DofFixedException::getValue
   @return Returns the value given in the constructor
*/

#endif
