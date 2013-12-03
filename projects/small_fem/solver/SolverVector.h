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

   Finaly, a SolverVector may be of the following scalar types:
   @li Real
   @li Complex
 */

template<typename scalar>
class SolverVector{
 private:
  size_t  N;
  scalar* v;

  mutable size_t      fail;
  mutable omp_lock_t* lock;

 public:
   SolverVector(void);
   SolverVector(size_t size);
  ~SolverVector(void);

  size_t getSize(void) const;
  size_t getNumberOfMutexWait(void) const;
  scalar get(size_t i) const;

  void    clear(void);
  void    resize(size_t size);
  void    add(size_t i, scalar value);
  scalar* getData(void);

 private:
  void init(size_t size);
};

/**
   @fn SolverVector::SolverVector()
   @param size An integer

   Instanciates a new SolverVector of size zero.
   **

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

   @fn SolverVector::get
   @param i A valid index of this SolverVector
   @return Returns a @em copy of the value of this SolverVector
   stored at the ith index (stating at index 0)
   **

   @fn SolverVector::clear
   Sets this SolverVector size to zero
   **

   @fn SolverVector::resize
   @param size An integer

   Clears this SolverVector and intialises a new one of the given size
   and with null values
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


//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SolverVectorInclusion.h"

//////////////////////
// Inline Functions //
//////////////////////

template<typename scalar>
inline
scalar SolverVector<scalar>::get(size_t i) const{
  // Temp
  scalar tmp;

  // Take mutex for the given index
  if(!omp_test_lock(&lock[i])){
    // If we can't we increase the total mutex wait ...
    fail++;

    // End we retry to take the lock
    omp_set_lock(&lock[i]);
  }

  // Add the given pair (column, value)
  tmp = v[i];

  // Free mutex
  omp_unset_lock(&lock[i]);

  // Return
  return tmp;
}

template<typename scalar>
inline
void SolverVector<scalar>::add(size_t i, scalar value){
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

#endif
