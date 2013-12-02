/////////////////////////////////////////////////
// Templates Implementations for SolverMatrix: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include <sstream>
#include <fstream>

template<typename scalar>
SolverMatrix<scalar>::SolverMatrix(void){
}

template<typename scalar>
SolverMatrix<scalar>::SolverMatrix(size_t nRow, size_t nCol){
  // Mutex wait //
  fail = 0;

  // Matrix size //
  this->nRow = nRow;
  this->nCol = nCol;

  // Matrix data //
  data = new std::list<std::pair<size_t, scalar> >[nRow];
  lock = new omp_lock_t[nRow];

  // Init locks //
  for(size_t i = 0; i < nRow; i++)
    omp_init_lock(&lock[i]);
}

template<typename scalar>
SolverMatrix<scalar>::~SolverMatrix(void){
  for(size_t i = 0; i < nRow; i++)
    omp_destroy_lock(&lock[i]);

  delete[] lock;
  delete[] data;
}

template<typename scalar>
size_t SolverMatrix<scalar>::nRows(void) const{
  return nRow;
}

template<typename scalar>
size_t SolverMatrix<scalar>::nColumns(void) const{
  return nCol;
}

template<typename scalar>
size_t SolverMatrix<scalar>::getNumberOfMutexWait(void) const{
  return fail;
}

template<typename scalar>
size_t SolverMatrix<scalar>::serialize(std::vector<int>&    rowVector,
                                       std::vector<int>&    colVector,
                                       std::vector<scalar>& valueVector){
  // Reduce the data vector such that we don't have redundant entries
  sortAndReduce();

  // Number of non zero entries
  size_t nNZ = 0;
  for(size_t i = 0; i < nRow; i++)
    nNZ += data[i].size();

  // Offset per row
  std::vector<size_t> offset(nRow);
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
    typename std::list<std::pair<size_t, scalar> >::iterator it;
    typename std::list<std::pair<size_t, scalar> >::iterator end;

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

template<typename scalar>
size_t SolverMatrix<scalar>::serializeCStyle(std::vector<int>&    rowVector,
                                             std::vector<int>&    colVector,
                                             std::vector<scalar>& valueVector){
  // Do it Fortran Style
  const size_t nNZ = serialize(rowVector, colVector, valueVector);

  // And substract 1 from indices
  for(size_t i = 0; i < nNZ; i++)
    rowVector[i]--;

  for(size_t i = 0; i < nNZ; i++)
    colVector[i]--;

  // Return
  return nNZ;
}

template<typename scalar>
std::string SolverMatrix<scalar>::toString(void) const{
  std::stringstream stream;
  typename std::list<std::pair<size_t, scalar> >::iterator it;
  typename std::list<std::pair<size_t, scalar> >::iterator end;

  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << "(" << i << ", " << it->first << "): " << it->second
             << std::endl;
  }

  return stream.str();
}

template<typename scalar>
void SolverMatrix<scalar>::writeToMatlabFile(std::string fileName,
                                             std::string matrixName) const{
  std::ofstream stream;
  stream.open(fileName.c_str());
  stream << toMatlab(matrixName) << std::endl;
  stream.close();
}

template<typename scalar>
void SolverMatrix<scalar>::sortAndReduce(void){
  // Each thread takes a set of rows
  #pragma omp parallel for
  for(size_t i = 0; i < nRow; i++){
    // Sort row[i]
    data[i].sort(sortPredicate);

    // Remove duplicate in row[i] by adding them
    typename std::list<std::pair<size_t, scalar> >::iterator start =
      data[i].begin();
    typename std::list<std::pair<size_t, scalar> >::iterator stop  =
      data[i].begin();
    typename std::list<std::pair<size_t, scalar> >::iterator end   =
      data[i].end();

    size_t idx;
    scalar val;

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
        typename std::list<std::pair<size_t, scalar> >::iterator tmp = stop;
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

template<typename scalar>
std::string SolverMatrix<scalar>::matlabCommon(std::string matrixName) const{
  // Stream
  typename std::list<std::pair<size_t, scalar> >::iterator it;
  typename std::list<std::pair<size_t, scalar> >::iterator end;
  std::stringstream stream;

  // Call to 'sparse'
  stream << matrixName << " = sparse(";

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

  // Return
  return stream.str();
}
