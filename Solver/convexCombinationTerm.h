// Gmsh - Copyright (C) 1997-2009 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _CONVEX_COMBINATION_TERM_H_
#define _CONVEX_COMBINATION_TERM_H_

#include <assert.h>
#include "femTerm.h"
#include "simpleFunction.h"
#include "Gmsh.h"
#include "GModel.h"
#include "SElement.h"
#include "fullMatrix.h"
#include "Numeric.h"

class convexCombinationTerm : public femTerm<double,double> {
 protected:
  const simpleFunction<double> *_k;
  const int _iField;
 public:
  convexCombinationTerm(GModel *gm, int iField, simpleFunction<double> *k)
    : femTerm<double,double>(gm), _iField(iField), _k(k) {}
  virtual int sizeOfR(SElement *se) const
  {
    return se->getMeshElement()->getNumVertices();
  }
  virtual int sizeOfC(SElement *se) const 
  { 
    return se->getMeshElement()->getNumVertices();
  }
  Dof getLocalDofR(SElement *se, int iRow) const
  {
    return Dof(se->getMeshElement()->getVertex(iRow)->getNum(), 
               Dof::createTypeWithTwoInts(0, _iField));
  }
  virtual void elementMatrix(SElement *se, fullMatrix<double> &m) const
  {

    MElement *e = se->getMeshElement();
    m.set_all(0.);
    const int nbNodes = e->getNumVertices();
    const double _diff = 1.0;
    for (int j = 0; j < nbNodes; j++){
      for (int k = 0; k < nbNodes; k++){
        m(j,k) = -1.*_diff;
      }
      m(j,j) = (nbNodes - 1) * _diff;
    }
  }


};

#endif
