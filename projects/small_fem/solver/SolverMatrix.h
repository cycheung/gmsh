#ifndef _SOLVERMATRIX_H_
#define _SOLVERMATRIX_H_

#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <omp.h>

/**
   @class SolverMatrix
   @brief A class handling a Solver linear system matrix

   This class represents a matrix used by a Solver.

   Once constructed, a SolverMatrix can be serialized into
   three vectors, r[], c[] and a[], such that a[k] is
   the entry of the sparse matrix at row (r[k] - 1) and column (c[k] - 1).

   During the construction of the matrix, multiple values
   can be added in a thread-safe manner.
*/

class SolverMatrix{
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
   SolverMatrix(size_t nRow, size_t nCol);
  ~SolverMatrix(void);

  size_t nRows(void) const;
  size_t nColumns(void) const;
  size_t getNumberOfMutexWait(void) const;

  void   add(size_t row, size_t col, double value);
  size_t serialize(std::vector<int>&    rowVector,
                   std::vector<int>&    colVector,
                   std::vector<double>& valueVector);

  size_t serializeCStyle(std::vector<int>&    rowVector,
                         std::vector<int>&    colVector,
                         std::vector<double>& valueVector);

  std::string toString(void) const;
  std::string toMatlab(std::string matrixName) const;
  void writeToMatlabFile(std::string fileName, std::string matrixName) const;

 private:
  SolverMatrix(void);

  void sortAndReduce(void);

  static bool sortPredicate(const std::pair<size_t, double>& a,
                            const std::pair<size_t, double>& b);
};

/**
   @fn SolverMatrix::SolverMatrix(size_t nRow, size_t nCol)
   @param nRow The number of row of this SolverMatrix
   @param nCol The number of column of this SolverMatrix

   Instanciates an new SolverMatrix of zero values with the given number of
   rows and columns
   **

   @fn SolverMatrix::~SolverMatrix

   Deletes this SolverMatrix
   **

   @fn SolverMatrix::nRows;
   @return Returns the number of rows of this SolverMatrix
   **

   @fn SolverMatrix::nColumns;
   @return Returns the number of columns of this SolverMatrix
   **

   @fn SolverMatrix::getNumberOfMutexWait
   @return Returns the total number of threads that were blocked by a mutex
   **

   @fn SolverMatrix::add
   @param row A row index of this SolverMatrix
   @param col A column index of this SolverMatrix
   @param value A real number

   Adds the given value at the given row and column
   **

   @fn SolverMatrix::serialize
   @param rowVector A vector of integers
   @param colVector A vector of integers
   @param valueVector A vector of reals

   @return Returns the number of non zero entries in this SolverMatrix

   Serializes this SolverMatrix.
   The three given vector will be such that:
   A[rowVector[k] - 1, colVector[k] - 1] = valueVector[k],
   where A is this SolverMatrix.
   **

   @fn SolverMatrix::serializeCStyle

   @param rowVector A vector of integers
   @param colVector A vector of integers
   @param valueVector A vector of reals

   @return Returns the number of non zero entries in this SolverMatrix

   Same as SolverMatrix::serialize, but does the job in C style indexing.
   Thus, the three given vector will be such that:
   A[rowVector[k], colVector[k]] = valueVector[k],
   **

   @fn SolverMatrix::toString
   @return Returns a string describing this SolverMatrix
   **

   @fn SolverMatrix::toMatlab
   @param matrixName A string
   @return Returns a string that can be used in Octave/Matlab
   to reproduce this SolverMatrix, whose name will be the given one
   **

   @fn SolverMatrix::writeToMatlabFile
   @param fileName A string
   @param matrixName A string

   Writes this matrix in Octave/Matlab format into the given file,
   and with the given name
 */

//////////////////////
// Inline Functions //
//////////////////////

inline size_t SolverMatrix::nRows(void) const{
  return nRow;
}

inline size_t SolverMatrix::nColumns(void) const{
  return nCol;
}

inline size_t SolverMatrix::getNumberOfMutexWait(void) const{
  return fail;
}

inline void SolverMatrix::add(size_t row, size_t col, double value){
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

inline bool SolverMatrix::sortPredicate(const std::pair<size_t, double>& a,
                                        const std::pair<size_t, double>& b){
  return a.first < b.first;
}


#endif
