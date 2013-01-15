#ifndef _WRITERMSH_H_
#define _WRITERMSH_H_

#include <fstream>

#include "BasisLagrange.h"
#include "Writer.h"

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

   The @em data of this Writer are either:
   @li A set of scalar (or vectorial) values, 
   defined on the @em vertices of the mesh, of @f$f@f$. 
   @li A System (or EigenSystem) solution

   A WriterMsh can write a .@c msh file of @f$f@f$

   @note
   If no data are given, WriteMsh will write a @c .msh file
   with @em only the given @em mesh.@n@n
   If no domain is given, an Exception will be thrown 
   when WriterMsh::write() is called.

   @todo
   Multi Topology for Hybrid Mesh (Adaptive View)@n
   Vectorial Adaptive View@n
*/

class WriterMsh: public Writer{
 protected:
  mutable std::ofstream* out;
  mutable BasisLagrange* lBasis;
  
 public:
  WriterMsh(void);
  
  virtual ~WriterMsh(void);
  
  virtual void write(const std::string name) const;

 protected:
  void writeHeader(void) const;
  void writeNodes(void) const;
  void writeElements(void) const;
  void writeInterpolationScheme(void) const;
  void writeNodalValuesHeader(const std::string name) const;
  void writeNodalValuesFromNode(void) const;
  void writeNodalValuesFromSol(void) const;  
  void writeNodalValuesFooter(void) const;
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
