#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <omp.h>

/**
   @class SparseMatrix
   @brief A class handling a sparse matrix

   This class represents a sparse matrix.

   Once constructed, a SparseMatrix can be serialized into
   three vectors, r[], c[] and a[], such that a[k] is
   the entry of the sparse matrix at row (r[k] - 1) and column (c[k] - 1).

   During the construction of the matrix, multiple values
   can be added in a thread-safe manner.
*/

class SparseMatrix{
 private:
  // Total Number of waiting mutex //
  size_t fail;

  // Size //
  size_t nRow;
  size_t nCol;

  // Data //
  // Each row is a list of the column entries
  std::list<std::pair<size_t, double> >* data;

  // Each row got also an omp lock to be thread-safe
  omp_lock_t* lock;

 public:
   SparseMatrix(size_t nRow, size_t nCol);
  ~SparseMatrix(void);

  size_t nRows(void) const;
  size_t nColumns(void) const;

  void   add(size_t row, size_t col, double value);
  size_t serialize(std::vector<int>&    rowVector,
                   std::vector<int>&    colVector,
                   std::vector<double>& valueVector);

  size_t      getNumberOfMutexWait(void) const;
  std::string toString(void) const;
  std::string toMatlab(void) const;

 private:
  SparseMatrix(void);

  void sortAndReduce(void);

  static bool sortPredicate(const std::pair<size_t, double>& a,
                            const std::pair<size_t, double>& b);
};

/**
   @fn SparseMatrix::SparseMatrix(size_t nRow, size_t nCol)
   @param nRow The number of row of this SparseMatrix
   @param nCol The number of column of this SparseMatrix

   Instanciates an new SparseMatrix of zero values with the given number of
   rows and columns
   **

   @fn SparseMatrix::~SparseMatrix

   Deletes this SparseMatrix
   **

   @fn SparseMatrix::nRows;
   @return Returns the number of rows of this SparseMatrix
   **

   @fn SparseMatrix::nColumns;
   @return Returns the number of columns of this SparseMatrix
   **

   @fn SparseMatrix::add
   @param row A row index of this SparseMatrix
   @param col A column index of this SparseMatrix
   @param value A real number

   Adds the given value at the given row and column
   **

   @fn SparseMatrix::serialize
   @param rowVector A vector of integers
   @param colVector A vector of integers
   @param valueVector A vector of reals

   @return Returns the number of non zero entries in this SparseMatrix

   Serializes this SparseMatrix.
   The three given vector will be such that:
   A[rowVector[k] - 1, colVector[k] - 1] = valueVector[k],
   where A is this SparseMatrix.
   **

   @fn SparseMatrix::getNumberOfMutexWait
   @return Returns the total number of threads that were blocked by a mutex
   **

   @fn SparseMatrix::toString
   @return Returns a string describing this SparseMatrix
   **

   @fn SparseMatrix::toMatlab
   @return Returns a string that can be used in Octave/Matlab
   to reproduce this SparseMatrix
   **
 */

/////////////////////
// Inline Function //
/////////////////////

inline size_t SparseMatrix::nRows(void) const{
  return nRow;
}

inline size_t SparseMatrix::nColumns(void) const{
  return nCol;
}

inline void SparseMatrix::add(size_t row, size_t col, double value){
  // Take mutex for the given row
  if(!omp_test_lock(&lock[row])){
    // If we can't we increase the total mutex wait ...
    fail++;

    // End we retry to take the lock
    omp_set_lock(&lock[row]);
  }

  // Add the given pair (column, value)
  data[row].push_back(std::pair<size_t, double>(col, value));

  // Free mutex
  omp_unset_lock(&lock[row]);
}

inline size_t SparseMatrix::getNumberOfMutexWait(void) const{
  return fail;
}

inline bool SparseMatrix::sortPredicate(const std::pair<size_t, double>& a,
                                        const std::pair<size_t, double>& b){
  return a.first < b.first;
}


#endif
