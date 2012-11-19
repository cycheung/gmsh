#ifndef _WRITER_H_
#define _WRITER_H_

#include "fullMatrix.h"
#include "MElement.h"
#include "MVertex.h"
#include "System.h"
#include "EigenSystem.h"

#include <string>
#include <vector>

/**
   @interface Writer
   @brief Common interface to Write Data into a file
   
   This is the common @em interface for Writer%s.@n

   A Writer is a class that can @em write a set of @em data
   into a @em file.@n

   Those data @em may be defined on a given @em domain.@n

   The @em exact meaning of the @em data and of the @em domain 
   @em must be specified by the actual @em implementation.

   @note
   A Writer is an @em interface, so it @em can't be instanciated
 */

class Writer{
 protected:
  bool ownSol;
  bool hasDomain;
  bool hasValue;
  bool isScalar;
  bool isNodal;

  int  N;
  int  E;

  const std::vector<const MElement*>* element;
  const std::vector<MVertex*>*        node;

  // Nodal Values //
  const std::vector<double>*              nodalScalarValue;
  const std::vector<fullVector<double> >* nodalVectorValue;

  // Values with inteprolation scheme //
  const FunctionSpace*      fs;
  const DofManager*         dofM;
  const fullVector<double>* sol;

 public:
  virtual ~Writer(void);

  virtual void write(const std::string name) const = 0;

  void setValues(const std::vector<double>& value);
  void setValues(const std::vector<fullVector<double> >& value);
  void setValues(const System& system);
  void setValues(const EigenSystem& system, unsigned int eigenNumber);

  void setDomain(const std::vector<const MElement*>& element);

 protected:
  Writer(void);

  static const fullVector<double>* 
    getSol(const std::vector<std::vector<std::complex<double> > >& eVector,
	   unsigned int eigenNumber);
};


/**
   @fn Writer::~Writer
   Deletes this Writer
   **

   @fn Writer::write
   @param name The name of the file to write into 
   (@em without extensions)
   
   Writes the Writer's Data into the given file
   **

   @fn void Writer::setValues(const std::vector<double>& value)
   @param value A set of value (double)
   
   Set this Writer's Data to the given values
   **

   @fn void Writer::setValues(const std::vector<fullVector<double> >& value)
   @param value A set of value (fullVector<double>)
   
   Set this Writer's Data to the given values
   **

   @fn void Writer::setValues(const System& value)
   @param value A System
   
   Use the System Solution for Data

   @warning
   Writer::setDomain() will be called with the Support
   of the FunctionSpace of the System
   **

   @fn Writer::setDomain
   @param element A set of MElement%s
   
   Set this Writer's Domain to the given elements
   **

   @internal
   @fn Writer::Writer
   This constructor 
   @endinternal
   **
 */

#endif
