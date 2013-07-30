#ifndef _WRITERMSH_H_
#define _WRITERMSH_H_

#include <fstream>

#include "BasisLagrange.h"
#include "Writer.h"

/**
   @class WriterMsh
   @brief A Writer for
   <a href="http://www.geuz.org/gmsh">gmsh</a> .msh file format.

   This class is a Writer implementing the
   <a href="http://www.geuz.org/gmsh">gmsh</a> .msh file format.

   Let @f$f@f$ be a scalar (or vectorial) function.

   The domain of this Writer is a set of MElement%s
   defining the geometrical support (a mesh) of @f$f@f$.

   The data of this Writer are either:
   @li A set of scalar (or vectorial) values
   @li A System (or EigenSystem) solution

   These data may be defined either:
   @li On the Vertices of each Element
   @li On each Element (Volume)

   A WriterMsh can write a .msh file of @f$f@f$

   If no data are given, WriteMsh will write a .msh file
   with only the given mesh.

   If no domain is given, an Exception will be thrown
   when WriterMsh::write() is called.
*/

class WriterMsh: public Writer{
 private:
  mutable std::ofstream* out;
  mutable BasisLagrange* lBasis;

 public:
  WriterMsh(void);

  virtual ~WriterMsh(void);

  virtual void write(const std::string name) const;
  virtual void write(const std::string name,
                     const std::string type) const;

 private:
  void writeInit(const std::string name) const;
  void writeFinalize(void) const;

  void writeAsNodalValues(const std::string name)  const;
  void writeAsVolumeValues(const std::string name) const;

  void writeHeader(void) const;
  void writeNodes(void) const;
  void writeElements(void) const;

  void writeInterpolationScheme(void) const;
  void writeNodalValuesHeader(const std::string name) const;
  void writeNodalValuesFromNode(void) const;
  void writeNodalValuesFromSol(void) const;
  void writeNodalValuesFooter(void) const;

  void writeVolumeValuesHeader(const std::string name) const;
  void writeVolumeValues(void) const;
  void writeVolumeValuesFooter(void) const;
};

/**
   @fn WriterMsh::WriterMsh
   Instantiates a new WriterMsh
   **

   @fn WriterMsh::~WriterMsh
   Deletes this WriterMsh
   **

   @fn WriterMsh::write(const std::string) const
   @param name The name of the file to write into (without extensions)

   Same as WriterMsh::write(name, vertex);
   **

   @fn WriterMsh::write(const std::string, const std::string) const
   @param name The name of the file to write into (without extensions)
   @param type A string

   Writes the Writer's Data into the given file

   If type is equal to:
   @li vertex, the given data will be used as nodal Values
   (that is defined on the vertices of an element)
   @li volume, the given data will be used as element Values
   (that is defined on an element)

   If no data are given, WriteMsh will write a .msh file
   with only the given mesh.

   If no domain is given, an Exception will be thrown.
*/

#endif
