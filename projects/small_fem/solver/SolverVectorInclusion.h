/////////////////////////////////////////////////
// Templates Implementations for SolverVector: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

template<typename scalar>
SolverVector<scalar>::SolverVector(void){
  // Alloc //
  fail = 0;
  N    = 0;
  v    = NULL;
  lock = NULL;
}

template<typename scalar>
SolverVector<scalar>::SolverVector(size_t size){
  init(size);
}

template<typename scalar>
SolverVector<scalar>::~SolverVector(void){
  clear();
}

template<typename scalar>
void SolverVector<scalar>::init(size_t size){
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
size_t SolverVector<scalar>::getSize(void) const{
  return N;
}

template<typename scalar>
size_t SolverVector<scalar>::getNumberOfMutexWait(void) const{
  return fail;
}

template<typename scalar>
void SolverVector<scalar>::clear(void){
  if(lock)
    delete[] lock;

  if(v)
    delete[] v;

  N    = 0;
  fail = 0;
  v    = NULL;
  lock = NULL;
}

template<typename scalar>
void SolverVector<scalar>::resize(size_t size){
  clear();
  init(size);
}

template<typename scalar>
scalar* SolverVector<scalar>::getData(void){
  return v;
}
