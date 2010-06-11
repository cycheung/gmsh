//
// C++ Interface: terms
//
// Description: Derivate class of SimpleFunction to include a time dependency
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMPLEFUNCTIONTIME_H_
#define SIMPLEFUNCTIONTIME_H_
template<class scalar>
class simpleFunctionTime : public simpleFunction<scalar>{
  protected :
    double time;
    scalar _val;
  public :
    simpleFunctionTime(scalar val=0,double t=1.) : simpleFunction<scalar>(val), _val(val), time(t){}; // time=1 by default to avoid set time for Static linear Scheme
    ~simpleFunctionTime(){};
    virtual scalar operator () (double x, double y, double z) const { return time*_val; }
    void setTime(const double t){time=t;}
};
#endif // SIMPLEFUNCTIONTIME
