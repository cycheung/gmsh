#ifndef _ELEMENTSOLUTION_H_
#define _ELEMENTSOLUTION_H_

#include <string>
#include <vector>

#include "GroupOfElement.h"
#include "PViewDataGModel.h"

/**
   @class ElementSolution
   @brief Handles per-element post-processing values

   This class associates to a group of element (GroupOfElement)
   a per-element value.

   A ElementSolution can write a post-processing map in the
   <a href="http://www.geuz.org/gmsh">gmsh</a> .msh file format.

   The same instance of ElementSolution may handle
   multiple sets of per-element values.
   Each set is associated to an integer called a 'step'.
   Each step is also associated to real value called 'time'.
   If two sets have the same step, the last one override to first one.
   Differents steps may have the same time.
 */

class ElementSolution{
 private:
  PViewDataGModel* pView;

 public:
   ElementSolution(void);
  ~ElementSolution(void);

  void clear(void);
  void addValues(size_t step,
                 double time,
                 const GroupOfElement& goe,
                 const std::vector<double>& value);

  void write(std::string fileName) const;
};


/**
   @fn ElementSolution::ElementSolution
   Instanciates a new ElementSolution which is empty
   (no elements and coefficients)
   **

   @fn ElementSolution::~ElementSolution
   Deletes this ElementSolution
   **

   @fn ElementSolution::clear
   This ElementSolution is now empty
   **

   @fn ElementSolution::addValues
   @param step An integer value
   @param time A real value
   @param goe A GroupOfElement
   @param value A set of values

   Adds the given set of values to this ElementSolution
   with the given step and time.

   We have that value(i) is the value of GroupOfElement::getAll()[i]
   **

   @fn ElementSolution::write
   @param fileName A file name (without extension)

   Writes this ElementSolution in
   <a href="http://www.geuz.org/gmsh">gmsh</a> .msh file format
   into the given file
 */

#endif
