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
#include "reduction.h"
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
 protected :
  int _num; // number of law (must be unique !)
 public :
  enum matname{linearElasticPlaneStress, linearElasticPlaneStressWithFracture};
  matname _type; // used to make a dynamic_cast to access to specific data of the law
  // constructor
  materialLaw(const int num) : _num(num){}
  ~materialLaw(){}
  materialLaw(const materialLaw &source){
    _num = source._num;
    _type = source._type;
  }
  materialLaw& operator=(const materialLaw &source){
    _num = source._num;
    _type = source._type;
    return *this;
  }
  virtual int getNum() const{return _num;}
  virtual matname getType() const{return _type;}
  virtual void setType(const matname t){_type=t;}
  static void registerBindings(binding *b){
    classBinding *cb = b->addClass<materialLaw>("materialLaw");
    cb->setDescription("base class for material Law");
  }
};

// class for linear elastic law
class linearElasticLawPlaneStress : public materialLaw{
  protected :
    const double _E; // YoungModulus //Store ??
    const double _nu; // Poisson ratio // Store ??
    const double C11;
    LinearElasticShellHookeTensor H;
  public :
    linearElasticLawPlaneStress(const int num, const double E, const double nu) : materialLaw(num),
                                                                                      _E(E), _nu(nu), C11(E/((1-nu)*(1+nu))){
      this->setType(materialLaw::linearElasticPlaneStress);
    }
    ~linearElasticLawPlaneStress(){}
    linearElasticLawPlaneStress( const linearElasticLawPlaneStress &source) : materialLaw(source), _E(source._E), _nu(source._nu), C11(source.C11), H(source.H){
    }
    virtual void stress(const LocalBasis *lb,const tab6 &eps,tab6 &sig){ //template scal vect tensor
        H.set(lb,C11,_nu);
        sig[0] = H(0,0,0,0)*eps[0]+(H(0,0,0,1)+H(0,0,1,0))*eps[3]+H(0,0,1,1)*eps[1];
        sig[1] = H(1,1,0,0)*eps[0]+(H(1,1,0,1)+H(1,1,1,0))*eps[3]+H(1,1,1,1)*eps[1];
        sig[3] = H(0,1,0,0)*eps[0]+(H(0,1,0,1)+H(0,1,1,0))*eps[3]+H(0,1,1,1)*eps[1];
        sig[2]=sig[4]=sig[5]=0.;
    }
    virtual double getYoung() const{return _E;}
    virtual double getPoisson() const{return _nu;}
    static void registerBindings(binding *b){
      classBinding *cb = b->addClass<linearElasticLawPlaneStress>("linearElasticLawPlaneStress");
      cb->setDescription("A linear elastic plane stress material law");
      methodBinding *cm;
      // Constructor
      cm = cb->setConstructor<linearElasticLawPlaneStress,int,double,double>();
      cm->setArgNames("num","E","nu",NULL);
      cm->setDescription("First arg is a unique law number. Second is the Young Modulus and third is the poisson coefficient");
    }
};

class linearElasticLawPlaneStressWithFracture : public linearElasticLawPlaneStress{
 protected :
  // const when LUA with 6 double
  double _Gc;
  double _sigmac;
  double _beta;
  double _mu;
 public :
  linearElasticLawPlaneStressWithFracture(const int num, const double E, const double nu, const double Gc, const double sigmac, const double beta, const double mu) :
                                          linearElasticLawPlaneStress(num,E,nu), _Gc(Gc), _sigmac(sigmac), _beta(beta), _mu(mu){
    this->setType(materialLaw::linearElasticPlaneStressWithFracture);
  }
  linearElasticLawPlaneStressWithFracture(const int num, const double E, const double nu) :
                                          linearElasticLawPlaneStress(num,E,nu), _Gc(0.), _sigmac(0.), _beta(0.), _mu(0.){
    this->setType(materialLaw::linearElasticPlaneStressWithFracture);
  }
  // set operation (needed because impossible to give 6 double with LUA why ?)
  void setGc(const double gc){_Gc=gc;}
  void setSigmac(const double sig){_sigmac = sig;}
  void setBeta(const double b){_beta =b;}
  void setMu(const double m){_mu = m;}

  // get operation
  double getGc() const{return _Gc;}
  double getSigmac() const{return _sigmac;}
  double getBeta() const{return _beta;}
  double getMu() const{return _mu;}
  void getCohesiveReduction(const reductionElement &m0, const reductionElement &n0, const double deltan,
                            const double deltan_max, const double deltat, const double deltat_max, const double deltac, const bool tension,
                            reductionElement &nhatmean, reductionElement &mhatmean) const{
    // for now Mxx and Nxx component (other = 0)
    // nhatmean and mhatmean are supposed to be initialized (change this when others components will be computed)
    nhatmean.setAll(0.);
    mhatmean.setAll(0.);
    // monotonic decreasing cohesive law (Camacho & Ortiz 1996)
    if(tension){ // tension case
      if(deltan<=deltac and deltan >=0.){
        double c;
        if( deltan >= deltan_max) // loading case
          c = 1.-deltan/deltac;
        else if(deltan > 0.) //unloading case
          c = deltan/deltan_max - deltan/deltac;
        mhatmean(1,1) = m0(1,1)*c; // Change this ??
        nhatmean(1,1) = n0(1,1)*c;
        mhatmean(0,1) = mhatmean(1,0) = m0(0,1)*c*sign(deltat);
        nhatmean(0,1) = nhatmean(1,0) = n0(0,1)*c*sign(deltat);
 //     printf("%lf %lf %lf\n",delta,mhatmean[1][1],M0);
      }
//      else if(deltan<0.)
//        Msg::Error("Deltan is <0 in a tension case !");
    }
    else{ // compression case
      if(deltat<=deltac){
        double c;
        if(fabs(deltat) >= fabs(deltat_max)) // loading case
          c = 1.-fabs(deltat)/deltac;
        else //unloading case
          c = fabs(deltat)/fabs(deltat_max) - fabs(deltat)/deltac;
        mhatmean(0,1) = mhatmean(1,0) = m0(0,1)*c*sign(deltat);
        nhatmean(0,1) = nhatmean(1,0) = n0(0,1)*c*sign(deltat);
      }
    }
  }
  static void registerBindings(binding *b){
    classBinding *cb = b->addClass<linearElasticLawPlaneStressWithFracture>("linearElasticLawPlaneStressWithFracture");
    cb->setDescription("A linear elastic plane stress material law with fracture possibility");
    methodBinding *cm;
    // Constructor
    cm = cb->setConstructor<linearElasticLawPlaneStressWithFracture,int,double,double>();
    cm->setArgNames("num","E","nu",NULL);
    cm->setDescription("Arguments are num E and nu. Others parameters must be pass thanks to set methods ");
    // methods
    cm = cb->addMethod("setGc", &linearElasticLawPlaneStressWithFracture::setGc);
    cm->setArgNames("gc",NULL);
    cm->setDescription("Set the value of fracture energy");
    cm = cb->addMethod("setSigmac", &linearElasticLawPlaneStressWithFracture::setSigmac);
    cm->setArgNames("sig",NULL);
    cm->setDescription("Set the value of fracture stress");
    cm = cb->addMethod("setBeta", &linearElasticLawPlaneStressWithFracture::setBeta);
    cm->setArgNames("b",NULL);
    cm->setDescription("Set the value of ratio : KII/KI");
    cm = cb->addMethod("setMu", &linearElasticLawPlaneStressWithFracture::setMu);
    cm->setArgNames("m",NULL);
    cm->setDescription("Set the value of friction coefficient");
  }
};

#endif //_MATERIALLAW_H_
