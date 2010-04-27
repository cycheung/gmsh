//
// Description : Parametric representation of a beam
//
//
// Author:   <Boris Sedji>,  04/2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef _BEAMPARAM_H_
#define _BEAMPARAM_H_

#include "SPoint2.h"
#include "MElement.h"

class SPoint2;

class BeamParam
{

  protected :

  public :

    virtual SPoint2 getPoint(double t) = 0;
    virtual BeamParam* cut(MElement *e) = 0;
    virtual double getSin() = 0;
    virtual double getCos() = 0;
    virtual double getC() = 0 ;
    virtual double getLength() = 0 ;
    virtual int getType() = 0;
    virtual std::vector <SPoint2>* getPoints() = 0;

};


#include "BeamParam.cpp"

#endif
