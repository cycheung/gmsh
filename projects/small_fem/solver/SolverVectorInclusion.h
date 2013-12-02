/////////////////////////////////////////////////
// Templates Implementations for SolverVector: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

template<typename scalar>
SolverVector<scalar>::SolverVector(void){
}

template<typename scalar>
SolverVector<scalar>::SolverVector(size_t size){
  // Alloc //
  fail = 0;
  N    = size;
  v    = new scalar[N];
  lock = new omp_lock_t[N];

  // Init //
  for(size_t i = 0; i < N; i++)
    v[i] = 0;

  for(size_t i = 0; i < N; i++)
    omp_init_lock(&lock[i]);

}

template<typename scalar>
SolverVector<scalar>::~SolverVector(void){
  delete[] lock;
  delete[] v;
}


template<typename scalar>
size_t SolverVector<scalar>::getSize(void) const{
  return N;
}

template<typename scalar>
size_t SolverVector<scalar>::getNumberOfMutexWait(void) const{
  return fail;
}

template<typename scalar>
scalar* SolverVector<scalar>::getData(void){
  return v;
}
