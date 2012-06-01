#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "Mesh.h"
#include "Formulation.h"
#include "fullMatrix.h"

#include "Interpolator.h"

#include <fstream>
#include <string>
#include <vector>

/**
   @class Solution
   @brief Writes the solution of a problem

   This class writes the solution of a finite element
   problem in a @c .msh 
   (<a href="http://www.geuz.org/gmsh">gmsh</a>
   file format)
   file.@n

   The problem is defined by a Mesh and a Formulation.@n
   The solution is the given by the Mesh Entity values.@n
   The Entity type to use is given by the Formulation.

   @todo
   May be use a general 'Problem' class ??@n
   Allow multiple fields
*/

class Solution{
 private:
  const Mesh* msh;
  const std::vector<Element*>* element;
  const std::vector<Node*>*    node;
  int N;
  int E;

  Interpolator* interp;

  std::vector<double>* nodalScalarValue;
  const std::vector<fullVector<double>*>* nodalVectorValue;

 public:
   Solution(const Mesh& mesh, const Formulation& formulation);
  ~Solution(void);

  void write(const std::string fileName,
	     const std::string name) const;

 private:
  void writeHeader(std::ofstream& out) const;
  void writeNodes(std::ofstream& out) const;
  void writeElements(std::ofstream& out) const;
  void writeNodalValues(std::ofstream& out, 
			bool isScalar,
			const std::string name) const;
};

/**
   @fn Solution::Solution
   @param mesh The Mesh to use
   @param formulation The Formulation to use
   @return Returns a new Solution with the given parameters

   @fn Solution::~Solution
   @return Deletes this Solution

   @fn Solution::write
   @param fileName The path to the @c .msh 
   to write the solution into
   @param name Name of the solution
   @return Writes the solution in the given @c .msh file
 */

#endif
