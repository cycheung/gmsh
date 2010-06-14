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
    bool timeDependency;
  public :
    simpleFunctionTime(scalar val=0,double t=1., bool td=true) : simpleFunction<scalar>(val), _val(val), time(t), timeDependency(td){}; // time=1 by default to avoid set time for Static linear Scheme
    ~simpleFunctionTime(){};
    virtual scalar operator () (double x, double y, double z) const { if(timeDependency) return time*_val; else return _val;}
    void setTime(const double t){time=t;}
};
#endif // SIMPLEFUNCTIONTIME
