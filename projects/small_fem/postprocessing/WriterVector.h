#ifndef _WRITERVECTOR_H_
#define _WRITERVECTOR_H_

#include <ostream>
#include "Writer.h"

/**
   @class WriterVector
   @brief A Writer for Raw Vector Data
  
   This class is a Writer for Raw Vector Data.

   These data can be @em scalar or @em vectorial.

   @note
   For this Writer, domain has no meaning.@n
   So, setting a domain has @em no @em effect.

   @note
   The special file name @em stdout will write the 
   vector into the standard outpout.
*/

class WriterVector: public Writer{
 public:
  WriterVector(void);
  
  virtual ~WriterVector(void);
  
  virtual void write(const std::string name) const;
  
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

   @fn void WriterVector::write(const std::string name) const
   @param name The name of the file to write into 
   (@em without extensions)
   
   Writes the Writer's Data into the given file

   @note
   The special file name @em stdout will write the 
   vector into the standard outpout.
   **   
*/

#endif
