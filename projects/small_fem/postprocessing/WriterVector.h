#ifndef _WRITERVECTOR_H_
#define _WRITERVECTOR_H_

#include <ostream>
#include "Writer.h"

/**
   @class WriterVector
   @brief A Writer for Raw Vector Data

   This class is a Writer for Raw Vector Data.

   These data can be scalar or vectorial.

   For this Writer, domain has no meaning.
   So, setting a domain has no effect.

   The special file name stdout will write the
   vector into the standard outpout.
*/

class WriterVector: public Writer{
 public:
  WriterVector(void);

  virtual ~WriterVector(void);

  virtual void write(const std::string name) const;
  virtual void write(const std::string name,
                     const std::string type) const;

 private:
  void write(std::ostream& stream) const;
};

/**
   @fn WriterVector::WriterVector
   Instantiates a new WriterVector
   **

   @fn WriterVector::~WriterVector
   Deletes this WriterVector
   **

   @fn void WriterVector::write(const std::string) const
   @param name The name of the file to write into (without extensions)

   Writes the Writer's Data into the given file

   The special file name stdout will write the vector into the standard outpout.
   **

   @fn void WriterVector::write(const std::string, const std::string) const
   @param name The name of the file to write into (without extensions)
   @param type A string

   Does the same as WriterVector::writer(const std::string),
   no matter the type field
*/

#endif
