#include "ThreadVector.h"

ThreadVector::ThreadVector(void){
}

ThreadVector::ThreadVector(size_t size){
  // Alloc //
  fail = 0;
  N    = size;
  v    = new double[N];
  lock = new omp_lock_t[N];

  // Init //
  for(size_t i = 0; i < N; i++)
    v[i] = 0;

  for(size_t i = 0; i < N; i++)
    omp_init_lock(&lock[i]);

}

ThreadVector::~ThreadVector(void){
  delete[] lock;
  delete[] v;
}
