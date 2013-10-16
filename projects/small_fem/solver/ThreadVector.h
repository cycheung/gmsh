#ifndef _THREADVECTOR_H_
#define _THREADVECTOR_H_

#include <cstring>
#include <omp.h>

/**
   @class ThreadVector
   @brief Thread-Safe Vector

   This class implements a thread-safe vector.

   @warning
   This shouldn't be used.
   I need to rework fullVector.
 */

class ThreadVector{
 private:
  size_t      fail;
  size_t      N;
  double*     v;
  omp_lock_t* lock;

 public:
   ThreadVector(size_t size);
  ~ThreadVector(void);

  size_t getSize(void) const;
  size_t getNumberOfMutexWait(void) const;

  void add(size_t i, double value);

  double* getData(void);

 private:
  ThreadVector(void);
};

/**
   @fn ThreadVector::ThreadVector(size_t)
   @param size An integer

   Instanciates a new ThreadVector of the given size and with null entries.
   **

   @fn ThreadVector::~ThreadVector

   Deletes this ThreadVector
   **

   @fn ThreadVector::getSize
   @return Returns this ThreadVector size
   **

   @fn ThreadVector::getNumberOfMutexWait
   @return Returns the total number of threads that were blocked by a mutex
   **

   @fn ThreadVector::add
   @param i An index of this ThreadVector
   @param value A real number

   Adds the given value at the given index of this ThreadVector
   **

   @fn ThreadVector::getData
   @return Returns a pointer to the memory segment holding the values
   of this ThreadVector
*/

//////////////////////
// Inline Functions //
//////////////////////

inline size_t ThreadVector::getSize(void) const{
  return N;
}

inline size_t ThreadVector::getNumberOfMutexWait(void) const{
  return fail;
}

inline void ThreadVector::add(size_t i, double value){
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

inline double* ThreadVector::getData(void){
  return v;
}

#endif
