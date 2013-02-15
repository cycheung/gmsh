#include "Term.h"

Term::Term(void){
  once = false;
}

Term::~Term(void){
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
