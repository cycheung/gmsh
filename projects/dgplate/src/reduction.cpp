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
#include "reduction.h"
inline double scaldot(const SVector3 &a, const SVector3 &b){
  double c=0.;
  for(int i=0;i<3;i++)
    c+= a[i]*b[i];
  return c;
}
static inline void diff(const SVector3 &a,const SVector3 &b,SVector3 &c){
  for(int i=0;i<3;i++)
    c[i]=a[i]-b[i];
}

double epsilongd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::GradType> &Grads, const int nbdof, const fullMatrix<double> &disp){
  int nbFF=Grads.size();
  double eps = 0.,u_g,u_d;
  SVector3 u;
  for(int i=0;i<nbFF;i++){
    u(0)=disp(i+nbdof,0); u(1)=disp(i+nbFF+nbdof,0); u(2)=disp(i+2*nbFF+nbdof,0);
    u_g=scaldot(u,lb->getphi0(gamma));
    u_d=scaldot(u,lb->getphi0(delta));
    eps+=(u_g*Grads[i](delta)+u_d*Grads[i](gamma));
  }
  eps/=2.;
  return eps;
}

double rhogd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::HessType> &hess, const int nbdof, const fullMatrix<double> &disp){
  int nbFF=hess.size();
  double rho = 0.,u_t;
  SVector3 u;
  for(int i=0;i<nbFF;i++){
    u(0)=disp(i+nbdof,0); u(1)=disp(i+nbFF+nbdof,0); u(2)=disp(i+2*nbFF+nbdof,0);
    u_t=scaldot(u,lb->gett0());
    rho-=(u_t*hess[i](gamma,delta));
  }
  return rho;
}

void stressReduction(const LinearElasticShellHookeTensor *H,const std::vector<TensorialTraits<double>::GradType> &Grads,const LocalBasis *lb,const fullMatrix<double> &disp, const int nbdof, std::vector<SVector3> &n){
  n[0](0)=0.;n[0](1)=0.;n[0](2)=0.;
  n[1](0)=0.;n[1](1)=0.;n[1](2)=0.;
  double eps_gd;
  for(int gamma=0;gamma<2;gamma++)
    for(int delta=0;delta<2;delta++)
    {
      eps_gd=epsilongd(gamma,delta,lb,Grads,nbdof,disp);
        for(int alpha=0;alpha<2;alpha++)
          for(int beta=0;beta<2;beta++)
            n[alpha](beta) += H->get(alpha,beta,gamma,delta)*eps_gd;
    }
}

void momentReduction(const LinearElasticShellHookeTensor *H,const std::vector<TensorialTraits<double>::HessType> &hess,const LocalBasis *lb,const fullMatrix<double> &disp, const int nbdof, std::vector<SVector3> &m){
  m[0](0)=0.;m[0](1)=0.;m[0](2)=0.;
  m[1](0)=0.;m[1](1)=0.;m[1](2)=0.;
  double rho_gd;
  for(int gamma=0;gamma<2;gamma++)
    for(int delta=0;delta<2;delta++)
    {
      rho_gd=rhogd(gamma,delta,lb,hess,nbdof,disp);
      for(int alpha=0;alpha<2;alpha++)
        for(int beta=0;beta<2;beta++)
          m[alpha](beta) += H->get(alpha,beta,gamma,delta)*rho_gd;
    }
}

void stressReductionHat(const std::vector<SVector3> &n,const LocalBasis *lb,std::vector<SVector3> &nhat){
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++){
      nhat[alpha](beta)=0.;
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++)
          nhat[alpha](beta)+=lb->getT(gamma,alpha)*n[gamma](delta)*lb->getT(delta,beta);
    }
}
// Should be somewhere else ??
void displacementjump(const std::vector<double> &Val_m,const int n_m,const std::vector<double> &Val_p,const int n_p,const fullMatrix<double> & disp,SVector3 &ujump){
  SVector3 up(0.,0.,0.),um(0.,0.,0.);
  for(int i=0;i<3;i++){
    for(int j=0;j<n_m;j++)
      um(i)+=Val_m[j]*disp(j+i*n_m,0);
    for(int j=0;j<n_p;j++)
      up(i)+=Val_p[j]*disp(j+i*n_p+3*n_m,0);
  }
  //diff(up,um,ujump);
  ujump = up-um;
};
