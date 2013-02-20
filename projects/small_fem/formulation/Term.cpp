#include "Term.h"

Term::Term(void){
  once = false;
}

Term::~Term(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;
}

double Term::getTermOutCache(unsigned int dofI,
                             unsigned int dofJ,
                             unsigned int elementId) const{
  unsigned int i   = 0;
  unsigned int ctr = elementId;
  unsigned int off = (*orientationStat)[0];

  for(; elementId >= off && i < nOrientation; i++){
    off += (*orientationStat)[i + 1];
    ctr -= (*orientationStat)[i];
  }

  once    = true;
  lastId  = elementId;
  lastI   = i;
  lastCtr = ctr;

  return (*aM[i])(ctr, dofI * nFunction + dofJ);
}

void Term::computeA(fullMatrix<double>**& bM,
                    fullMatrix<double>**& cM){
  // Alloc //
  aM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    aM[s] = new fullMatrix<double>((*orientationStat)[s], nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++)
    // GEMM doesn't like matrices with 0 Elements
    if((*orientationStat)[s])
      aM[s]->gemm(*bM[s], *cM[s]);
}

void Term::clean(fullMatrix<double>**& bM,
                 fullMatrix<double>**& cM){

  for(unsigned int s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;
}
