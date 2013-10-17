#include "Term.h"

Term::Term(void){
  // One cache per thread
  const size_t nThread = omp_get_max_threads();

  once    = new bool[nThread];
  lastId  = new size_t[nThread];
  lastI   = new size_t[nThread];
  lastCtr = new size_t[nThread];

  // We did not get into the cache
  for(size_t i = 0; i < nThread; i++)
    once[i] = false;
}

Term::~Term(void){
  delete[] once;
  delete[] lastId;
  delete[] lastI;
  delete[] lastCtr;

  for(size_t s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;
}

double Term::getTermOutCache(size_t dofI,
                             size_t dofJ,
                             size_t elementId,
                             size_t threadId) const{
  size_t i   = 0;
  size_t ctr = elementId;
  size_t off = (*orientationStat)[0];

  for(; elementId >= off && i < nOrientation; i++){
    off += (*orientationStat)[i + 1];
    ctr -= (*orientationStat)[i];
  }

  once[threadId]    = true;
  lastId[threadId]  = elementId;
  lastI[threadId]   = i;
  lastCtr[threadId] = ctr;

  return (*aM[i])(ctr, dofI * nFunction + dofJ);
}

void Term::allocA(size_t nFunction){
  // Alloc //
  aM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    aM[s] = new fullMatrix<double>((*orientationStat)[s], nFunction);
}

void Term::computeA(fullMatrix<double>**& bM,
                    fullMatrix<double>**& cM){
  // Fill //
  for(size_t s = 0; s < nOrientation; s++)
    // GEMM doesn't like matrices with 0 Elements
    if((*orientationStat)[s])
      aM[s]->gemm(*bM[s], *cM[s]);
}

void Term::clean(fullMatrix<double>**& bM,
                 fullMatrix<double>**& cM){

  for(size_t s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(size_t s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;
}
