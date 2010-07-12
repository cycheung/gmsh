//
// C++ Interface: terms
//
// Description: Function needed to compute element of reduction n^alpha and m^alpha for thin bodies
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef REDUCTION_H_
#define REDUCTION_H_
#include "LinearElasticShellHookeTensor.h"
#include "LocalBasis.h"
#include "SVector3.h"

// class for store reduction element m(alpha,beta) and n(alpha,beta)
class reductionElement{
 protected :
   double mat[2][2];
 public :

  // constructors and destructor
  reductionElement(){for(int i=0;i<2;i++) for(int j=0;j<2;j++) mat[i][j]=0.;}
  ~reductionElement(){}
  reductionElement(const reductionElement &source){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
        mat[i][j] = source.mat[i][j];
  }
  reductionElement & operator = (const reductionElement &source){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
        mat[i][j] = source.mat[i][j];
    return *this;
  }
  // access to data
  double& operator() (const int i, const int j){return mat[i][j];}
  double operator() (const int i, const int j) const {return mat[i][j];}
  void setAll(const double d){for(int i=0;i<2;i++) for(int j=0;j<2;j++) mat[i][j]=0.;}
};

// class for stressTensor
class stressTensor;
static void matvectprod(const stressTensor &st, const SVector3 &v1, SVector3 &v2);
class stressTensor : public STensor3{
 public :
  stressTensor(const double v =0.0) : STensor3(v){}
  stressTensor(const reductionElement &red, const LocalBasis *lbs, const double v =0.0) : STensor3(v){
    STensor3 temp;
    for(int alpha = 0; alpha<2;alpha++)
      for(int beta=0; beta<2;beta++){
        tensprod(lbs->getphi0(alpha),lbs->getphi0(beta),temp);
        temp *= red(alpha,beta);
        *this += temp;
      }
  }
  ~stressTensor(){}
  stressTensor(const stressTensor & source) : STensor3(source){}
  stressTensor &operator= (const stressTensor &source){
    this->STensor3::operator=(source);
    return *this;
  }
  double getComponent(const SVector3 &n1, const SVector3 & n2) const{
    SVector3 temp(0., 0., 0.);
    matvectprod(*this,n2,temp);
    return dot(temp,n1);
  }
};

static void matvectprod(const stressTensor &st, const SVector3 &v1, SVector3 &v2){
  SVector3 temp;
  for(int i=0;i<3;i++){
    temp(0) = st(i,0); temp(1) = st(i,1); temp(2) = st(i,2);
    v2(i) = dot(temp,v1);
  }
}

static inline double scaldot(const SVector3 &a, const SVector3 &b);

static inline void diff(const SVector3 &a,const SVector3 &b,SVector3 &c);

double epsilongd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::GradType> &Grads,
                  const std::vector<double> &disp);

double rhogd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::GradType> &Grads,
             const std::vector<TensorialTraits<double>::HessType> &hess,const std::vector<double> &disp);

void stressReductionHat(const reductionElement &n,const LocalBasis *lb, reductionElement &nhat);
// Should be somewhere else ??
void displacementjump(const std::vector<double> &Val_m,const int n_m,const std::vector<double> &Val_p,const int n_p,
                      const std::vector<double> & disp,SVector3 &ujump);
void displacementjump(const std::vector<double> &Val_m,const int n_m,const std::vector<double> &Val_p,
                      const int n_p,const std::vector<double> & dispm,const std::vector<double> & dispp,SVector3 &ujump);
void rotationjump(const LocalBasis *lbs,const std::vector<SVector3> &Gradm,const int nbFFm, const int nbdofm, const LocalBasis *lbm,
                  const std::vector<SVector3> &Gradp,const int nbFFp,const LocalBasis *lbp,
                  const std::vector<double> &disp, double rjump[3]);
void rotationjump(const LocalBasis *lbs,const std::vector<SVector3> &Gradm,const int nbFFm, const LocalBasis *lbm,
                  const std::vector<SVector3> &Gradp,const int nbFFp,const LocalBasis *lbp,
                  const std::vector<double> &dispm, const std::vector<double> &dispp, double rjump[3]);
 #endif // reduction.h
