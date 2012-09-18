#ifndef _EXACTSOLUTION_H_
#define _EXACTSOLUTION_H_

#include <string>
#include <vector>

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "Writer.h"

/**
   @class ExactSolution
   @brief Exact Solution of a Problem

   This class can compute the @em exact solution of a given Problem.@n
   
   The problem is formulated by a @em derived class.

   @note
   This class can't be instanciated
*/


class ExactSolution{
 protected:
  bool isScalar;
  
  const GroupOfElement* domain;

  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;

 public:
  virtual ~ExactSolution(void);
  
  void write(const std::string name, Writer& writer) const;

 protected:
  ExactSolution(void);
  virtual void compute(void);

  virtual double             fScalar(double x, double y, double z);
  virtual fullVector<double> fVector(double x, double y, double z);
};


/**
   @internal
   @fn ExactSolution::ExactSolution
   @param goe The Meshed Domain

   Instanciates a new ExactSolution
   @endinternal
   **

   @fn ExactSolution::~ExactSolution
   Deletes this ExactSolution   
   **

   @fn ExactSolution::write
   @param name The file (@em without extension) where the 
   solution will be written
   @param writer The Writer to use

   Writes the solution into the given file@n

   The file format is given by the Writer
   **
*/

#endif
