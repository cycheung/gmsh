//
// Description : Parametric representation of a beam
//
//
// Author:   <Boris Sedji>,  04/2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef _BEAMTERM_H_
#define _BEAMTERM_H_

#include "BeamParam.h"
#include "MElement.h"
#include "Gauss.h"
#include "fullMatrix.h"
#include "functionSpace.h"

class BeamTerm
{

  protected :

  BeamParam *_Beam;
  MElement *_e;
  FunctionSpace<SVector3> *_LagSpace;

  public :

  BeamTerm(FunctionSpace<SVector3> *LagSpace,MElement *e, BeamParam *Beam);
  //void get(int npts,IntPt *GP,fullMatrix<double> &m);
  void get(fullMatrix<double> &m);
};


#include "BeamTerm.cpp"

#endif

