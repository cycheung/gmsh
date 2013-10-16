#include "sstream"
#include "SparseMatrix.h"

using namespace std;

SparseMatrix::SparseMatrix(void){
}

SparseMatrix::SparseMatrix(size_t nRow, size_t nCol){
  // Mutex wait //
  fail = 0;

  // Matrix size //
  this->nRow = nRow;
  this->nCol = nCol;

  // Matrix data //
  data = new list<pair<size_t, double> >[nRow];
  lock = new omp_lock_t[nRow];

  // Init locks //
  for(size_t i = 0; i < nRow; i++)
    omp_init_lock(&lock[i]);
}

SparseMatrix::~SparseMatrix(void){
  for(size_t i = 0; i < nRow; i++)
    omp_destroy_lock(&lock[i]);

  delete[] lock;
  delete[] data;
}

size_t SparseMatrix::serialize(vector<int>&    rowVector,
                               vector<int>&    colVector,
                               vector<double>& valueVector){
  // Reduce the data vector such that we don't have redundant entries
  sortAndReduce();

  // Number of non zero entries
  size_t nNZ = 0;
  for(size_t i = 0; i < nRow; i++)
    nNZ += data[i].size();

  // Offset per row
  vector<size_t> offset(nRow);
  offset[0] = 0;

  for(size_t i = 1, j = 0; i < nRow; i++, j++)
    offset[i] = offset[j] + data[j].size();

  // Allocate
  rowVector.resize(nNZ);
  colVector.resize(nNZ);
  valueVector.resize(nNZ);

  // Launch parallel portion
  #pragma omp parallel
  {
    // Temp Private Data
    list<pair<size_t, double> >::iterator it;
    list<pair<size_t, double> >::iterator end;

    // Row Vector
    #pragma omp for
    for(size_t i = 0; i < nRow; i++){
      it  = data[i].begin();
      end = data[i].end();

      for(size_t j = 0; it != end; it++, j++)
        rowVector[offset[i] + j] = (int)(i) + 1;
    }

    // Column Vector
    #pragma omp for
    for(size_t i = 0; i < nRow; i++){
      it  = data[i].begin();
      end = data[i].end();

      for(size_t j = 0; it != end; it++, j++)
        colVector[offset[i] + j] = (int)(it->first) + 1;
    }

    // Value Vector
    #pragma omp for
    for(size_t i = 0; i < nRow; i++){
      it  = data[i].begin();
      end = data[i].end();

      for(size_t j = 0; it != end; it++, j++)
        valueVector[offset[i] + j] = it->second;
    }
  }

  // Return number of non zero entries
  return nNZ;
}

string SparseMatrix::toString(void) const{
  stringstream stream;
  list<pair<size_t, double> >::iterator it;
  list<pair<size_t, double> >::iterator end;

  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << "(" << i << ", " << it->first << "): " << it->second << endl;
  }

  return stream.str();
}

string SparseMatrix::toMatlab(void) const{
  // Init
  stringstream stream;
  list<pair<size_t, double> >::iterator it;
  list<pair<size_t, double> >::iterator end;

  // Call to 'sparse'
  stream << "A = sparse(";

  // Rows
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << i + 1 << ", ";
  }
  stream << "], ";

  // Columns
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << it->first + 1 << ", ";
  }
  stream << "], ";

  // Values
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << std::scientific << it->second << ", ";
  }
  stream << "], ";

  // Number of rows and columns
  stream << nRow << ", " << nCol << ")";

  // Return
  return stream.str();
}

void SparseMatrix::sortAndReduce(void){
  // Each thread takes a set of rows
  #pragma omp parallel for
  for(size_t i = 0; i < nRow; i++){
    // Sort row[i]
    data[i].sort(sortPredicate);

    // Remove duplicate in row[i] by adding them
    list<pair<size_t, double> >::iterator start = data[i].begin();
    list<pair<size_t, double> >::iterator stop  = data[i].begin();
    list<pair<size_t, double> >::iterator end   = data[i].end();

    size_t idx;
    double val;

    if(start != end){
      idx = start->first;
      val = start->second;
      stop++;
    } // Not needed but compiler cannot see that it is redundant with while loop

    while(start != end){
      // If start == stop, we need to resart
      if(start == stop){
        idx = start->first;
        val = start->second;
        stop++;
      }

      // If we get a element with the same column index:
      // we continue to accumate in 'val'
      if(stop != end && stop->first == idx){
        val += stop->second;
        stop++;
      }

      // Else, we earse the superflous element
      else{
        // Create an iterator to the previous element
        list<pair<size_t, double> >::iterator tmp = stop;
        tmp--;

        if(tmp != start){
          // We need to erase superflous elements
          // But first, set 'tmp' to value accumulator
          tmp->second = val;

          // Then erase from start to tmp (excluded)
          data[i].erase(start, tmp);
        }

        // Restart from 'stop'
        start = stop;
      }
    }
  }
}
