#ifndef _WRITER_H_
#define _WRITER_H_

#include "fullMatrix.h"
#include "GroupOfElement.h"
#include "MElement.h"
#include "MVertex.h"
#include "System.h"
#include "SystemEigen.h"

#include <string>
#include <vector>

/**
   @interface Writer
   @brief Common interface to Write Data into a file

   This is the common interface for Writer%s.

   A Writer is a class that can write a set of data
   into a file.

   Those data may be defined on a given domain.

   The exact meaning of the data and of the domain
   must be specified by the actual implementation.
 */

class Writer{
 protected:
  bool ownSol;
  bool hasDomain;
  bool hasValue;
  bool isScalar;
  bool isNodal;

  size_t N;
  size_t E;

  const std::vector<const MElement*>* element;
  const std::vector<MVertex*>*        node;

  // Nodal Values //
  const std::vector<double>*              scalarValue;
  const std::vector<fullVector<double> >* vectorValue;

  // Values with inteprolation scheme //
  const FunctionSpace*      fs;
  const DofManager*         dofM;
  const fullVector<double>* sol;

 public:
  virtual ~Writer(void);

  virtual void write(const std::string name) const = 0;
  virtual void write(const std::string name,
                     const std::string type) const = 0;

  void setValues(const std::vector<double>& value);
  void setValues(const std::vector<fullVector<double> >& value);
  void setValues(const System& system);
  void setValues(const SystemEigen& system, size_t eigenNumber);

  void setDomain(const GroupOfElement& domain);

 protected:
  Writer(void);

  static const fullVector<double>*
    getSol(const std::vector<fullVector<std::complex<double> > >& eVector,
           size_t eigenNumber);
};


/**
   @fn Writer::~Writer
   Deletes this Writer
   **

   @fn Writer::write(const std::string) const
   @param name The name of the file to write into (without extensions)

   Writes the Writer's Data into the given file

   If the Data may be interpreted, by a particular Writer,
   in multiple way, the method uses the default choice
   **

   @fn Writer::write(const std::string, const std::string) const
   @param name The name of the file to write into (without extensions)
   @param type A string

   Writes the Writer's Data into the given file

   If the Data may be interpreted, by a particular Writer,
   in multiple way, the filed type is used to overcome this situation
   **

   @fn void Writer::setValues(const std::vector<double>& value)
   @param value A set of value (double)

   Sets this Writer's Data to the given value
   **

   @fn void Writer::setValues(const std::vector<fullVector<double> >& value)
   @param value A set of value (fullVector<double>)

   Sets this Writer's Data to the given value
   **

   @fn void Writer::setValues(const System&)
   @param system A System

   Uses the solution of the given System as Data

   Writer::setDomain() will be called with FunctionSpace::getSupport()
   as argument. The FunctionSpace comes from SystemAbstract::getFunctionSpace()
   **

   @fn void Writer::setValues(const SystemEigen& system, size_t)
   @param system An Eigen System (SystemEigen)
   @param eigenNumber A natural number

   Uses an Eigenvector of the given Eigen System as Data.
   The Eigenvector is defined by eigenNumber.

   Writer::setDomain() will be called with FunctionSpace::getSupport()
   as argument. The FunctionSpace comes from SystemAbstract::getFunctionSpace()
   **

   @fn Writer::setDomain
   @param domain A GroupOfElement

   Set this Writer's Domain with the given group of elements (GroupOfElement)
   **

   @internal
   @fn Writer::Writer
   This constructor
   @endinternal
   **
 */

#endif
