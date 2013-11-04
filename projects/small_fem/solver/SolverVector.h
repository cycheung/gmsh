#ifndef _SOLVERVECTOR_H_
#define _SOLVERVECTOR_H_

#include <cstring>
#include <omp.h>

/**
   @class SolverVector
   @brief A class handling a Solver linear system vector

   This class represents a vector used by a Solver.

   During the construction of the vector, multiple values
   can be added in a thread-safe manner.
 */

class SolverVector{
 private:
  size_t      fail;
  size_t      N;
  double*     v;
  omp_lock_t* lock;

 public:
   SolverVector(size_t size);
  ~SolverVector(void);

  size_t getSize(void) const;
  size_t getNumberOfMutexWait(void) const;

  void add(size_t i, double value);

  double* getData(void);

 private:
  SolverVector(void);
};

/**
   @fn SolverVector::SolverVector(size_t)
   @param size An integer

   Instanciates a new SolverVector of the given size and with null entries.
   **

   @fn SolverVector::~SolverVector

   Deletes this SolverVector
   **

   @fn SolverVector::getSize
   @return Returns this SolverVector size
   **

   @fn SolverVector::getNumberOfMutexWait
   @return Returns the total number of threads that were blocked by a mutex
   **

   @fn SolverVector::add
   @param i An index of this SolverVector
   @param value A real number

   Adds the given value at the given index of this SolverVector
   **

   @fn SolverVector::getData
   @return Returns a pointer to the memory segment holding the values
   of this SolverVector
*/

//////////////////////
// Inline Functions //
//////////////////////

inline size_t SolverVector::getSize(void) const{
  return N;
}

inline size_t SolverVector::getNumberOfMutexWait(void) const{
  return fail;
}

inline void SolverVector::add(size_t i, double value){
  // Take mutex for the given index
  if(!omp_test_lock(&lock[i])){
    // If we can't we increase the total mutex wait ...
    fail++;

    // End we retry to take the lock
    omp_set_lock(&lock[i]);
  }

  // Add the given pair (column, value)
  v[i] += value;

  // Free mutex
  omp_unset_lock(&lock[i]);
}

inline double* SolverVector::getData(void){
  return v;
}

#endif
