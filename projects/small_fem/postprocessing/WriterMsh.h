#ifndef _WRITERMSH_H_
#define _WRITERMSH_H_

#include "Writer.h"
#include <fstream>

/**
   @class WriterMsh
   @brief A Writer for 
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   .@c msh file format.
  
   This class is a Writer implementing the    
   <a href="http://www.geuz.org/gmsh">gmsh</a>
   .@c msh file format.

   Let @f$f@f$ be a scalar (or vectorial) function.

   The @em domain of this Writer is a set of MElement%s
   defining the geometrical support (a mesh) of @f$f@f$. 

   The @em data of this Writer are the set of scalar (or vectorial)
   values, defined on the @em vertices of the mesh, of @f$f@f$. 

   A WriterMsh can write a .@c msh file of @f$f@f$

   @note
   If no data are given, WriteMsh will write a @c .msh file
   with @em only the given @em mesh.@n@n
   If no domain is given, an Exception will be thrown 
   when WriterMsh::write() is called.
*/

class WriterMsh: public Writer{
 protected:
  mutable std::ofstream* out;
  
 public:
  WriterMsh(void);
  
  virtual ~WriterMsh(void);
  
  virtual void write(const std::string name) const;

 protected:
  void writeHeader(void) const;
  void writeNodes(void) const;
  void writeElements(void) const;
  void writeNodalValues(const std::string name) const;
};

/**
   @fn WriterMsh::WriterMsh
   Instantiates a new WriterMsh
   **

   @fn WriterMsh::~WriterMsh
   Deletes this WriterMsh
   **

   @fn WriterMsh::write
   @param name The name of the file to write into 
   (@em without extensions)
   
   Writes the Writer's Data into the given file

   @note
   If no data are given, WriteMsh will write a @c .msh file
   with @em only the given @em mesh.@n@n
   If no domain is given, an Exception will be thrown 
   when WriterMsh::write() is called.
   **   
*/

#endif
