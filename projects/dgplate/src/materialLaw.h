//
// C++ Interface: terms
//
// Description: Define material law
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _MATERIALLAW_H_
#define _MATERIALLAW_H_
#include "LocalBasis.h"
#include "LinearElasticShellHookeTensor.h"
// class with all material laws

// Problem when it puts in IPState.h ????
// type array for store in vector
struct tab6{
  protected :
    double a[6];
  public :
    tab6(){for(int i=0;i<6;i++) a[i]=0.;}
    ~tab6(){};
    tab6& operator= (tab6 source) {for(int i=0;i<6;i++) a[i]=source[i];}
    double& operator[](const int i){return a[i];}
    double operator[](const int i)const {return a[i];}
};


class materialLaw{
  public :
    enum matname{linearElasticPlaneStress, linearElasticPlaneStressWithFracture};
    virtual void bidon(){}; // One virtual function to use polymorphism
    //virtual void stress(const LocalBasis*, const double[6],double[6])=0;
};

// class for linear elastic law
class linearElasticLawPlaneStress : public materialLaw{
  protected :
    const double _E; // YoungModulus //Store ??
    const double _nu; // Poisson ratio // Store ??
    const double C11;
    LinearElasticShellHookeTensor H;
  public :
    linearElasticLawPlaneStress(const double E, const double nu) : _E(E), _nu(nu), C11(E/((1-nu)*(1+nu))){
    }
    virtual void stress(const LocalBasis *lb,const tab6 &eps,tab6 &sig){ //template scal vect tensor
        H.set(lb,C11,_nu);
        sig[0] = H(0,0,0,0)*eps[0]+(H(0,0,0,1)+H(0,0,1,0))*eps[3]+H(0,0,1,1)*eps[1];
        sig[1] = H(1,1,0,0)*eps[0]+(H(1,1,0,1)+H(1,1,1,0))*eps[3]+H(1,1,1,1)*eps[1];
        sig[3] = H(0,1,0,0)*eps[0]+(H(0,1,0,1)+H(0,1,1,0))*eps[3]+H(0,1,1,1)*eps[1];
        sig[2]=sig[4]=sig[5]=0.;
    }
};

class linearElasticLawPlaneStressWithFracture : public linearElasticLawPlaneStress{
 protected :
  const double _Gc;
  const double _sigmac;
 public :
  linearElasticLawPlaneStressWithFracture(const double E, const double nu, const double Gc, const double sigmac) :
                                          linearElasticLawPlaneStress(E,nu), _Gc(Gc), _sigmac(sigmac){}
  // get operation
  double getGc() const{return _Gc;}
  double getSigmac() const{return _sigmac;}
  void getCohesiveReduction(const double M0, const double N0, const double delta, const double delta_max,
                              const double deltac,std::vector<SVector3> &nhatmean, std::vector<SVector3> &mhatmean) const{
    // for now Mxx and Nxx component (other = 0)
    // nhatmean and mhatmean are supposed to be initialized (change this when others components will be computed)
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++){
        mhatmean[i][j]=0.;
        nhatmean[i][j]=0.;
      }
    // monotonic decreasing cohesive law
    if((0.<=delta) and (delta<=deltac)){
      double c;
      if(delta >= delta_max) // loading case
        c = 1.-delta/deltac;
      else //unloading case
        c = delta/delta_max - delta/deltac;
      mhatmean[1][1] = M0*c; // Change this ??
      nhatmean[0][0] = N0*c;
 //     printf("%lf %lf %lf\n",delta,mhatmean[1][1],M0);
    }
//    else if(delta> deltac)
//      printf("full broken\n");
//    else {printf("haha\n");Msg::Error("Case delta < 0 is not yet implemented for linearElasticLawPlaneStressWithFracture\n");}
  }
};

#endif //_MATERIALLAW_H_
