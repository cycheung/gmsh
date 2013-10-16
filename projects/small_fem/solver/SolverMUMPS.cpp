#include "SolverMUMPS.h"
#include "Exception.h"
#include "dmumps_c.h"

using namespace std;

SolverMUMPS::SolverMUMPS(void){
}

SolverMUMPS::~SolverMUMPS(void){
}

void SolverMUMPS::solve(SparseMatrix& A,
                        fullVector<double>& rhs,
                        fullVector<double>& x){
  // Is the given matrix square ? //
  const int size = A.nRows();

  if((size_t)(size) != A.nColumns())
    throw Exception("SolverMUMPS -- The given matrix is not square: (%d, %d)",
                    size, A.nColumns());

  // Is the right hand side of the good size ? //
  if(size != rhs.size())
    throw Exception("%s -- %s: %d (matrix is %d)",
                    "SolverMUMPS", "The given rhs do not have the right size",
                    rhs.size(), size);

  // Serialize the given matrix //
  vector<int>    row;
  vector<int>    col;
  vector<double> value;
  const int nNZ = A.serialize(row, col, value);

  // Init MUMPS //
  DMUMPS_STRUC_C id;

  id.job          =      -1; // Initialize MUMPS instance
  id.par          =       1; // Host processor participates to the job
  id.sym          =       0; // Unsymmetric matrix
  id.comm_fortran = -987654; // Use MPI COMM WORLD

  dmumps_c(&id);             // Do what is told in struct 'id'

  // Define the problem //
  id.icntl[4]  = 0;      // Matrix in assembled format
  id.icntl[17] = 0;      // Matrix is centralized on the host

  id.n   = size;         // Size of the (square) matrix of unknown
  id.nz  = nNZ;          // Number of non zero entries in the matrix
  id.irn = row.data();   // Row vector
  id.jcn = col.data();   // Column vector
  id.a   = value.data(); // Value vector
  id.rhs = &rhs(0);      // Right hand side

  // Output Settings //
  id.icntl[0] = -1;  // No Output
  id.icntl[1] = -1;  // ---------
  id.icntl[2] = -1;  // ---------
  id.icntl[3] = -1;  // ---------

  // Call MUMPS //
  id.job = 6; // Do the default steps to solve
  dmumps_c(&id);

  // The Right hand side is now the solution //
  // Turn the 'x' vector into a proxy of 'rhs'
  x.setAsProxy(rhs, 0, rhs.size());

  // Terminate MUMPS instance //
  id.job = -2; // Destroy MUMPS Instance
  dmumps_c(&id);
}
