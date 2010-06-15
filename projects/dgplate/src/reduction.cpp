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
#include <math.h>
// Compute the normal variation
void normalVariation(const LocalBasis *lb, const int nbFF, const std::vector<TensorialTraits<double>::GradType> &Grads,
                     const std::vector<double> &disp, SVector3 &dt, const int nbdofm=0){
  for(int i=0;i<3;i++) dt[i]=0.;
  SVector3 v1(0.,0.,0.);
  SVector3 v2(0.,0.,0.);
  v1 = crossprod(lb->gett0(),lb->getphi0(0));
  v2 = crossprod(lb->gett0(),lb->getphi0(1));

  SVector3 du1(0.,0.,0.);
  SVector3 du2(0.,0.,0.);
  for(int i=0;i<nbFF;i++){
    du1[0] += Grads[i](0) * disp[i+nbdofm]; du1[1] += Grads[i](0) * disp[i+nbFF+nbdofm]; du1[2] += Grads[i](0) * disp[i+nbFF+nbFF+nbdofm];
    du2[0] += Grads[i](1) * disp[i+nbdofm]; du2[1] += Grads[i](1) * disp[i+nbFF+nbdofm]; du2[2] += Grads[i](1) * disp[i+nbFF+nbFF+nbdofm];
  }
  dt = crossprod(lb->getphi0(0),du2) + lb->gett0() * scaldot(du1,v2) - ( crossprod(lb->getphi0(1),du1) + lb->gett0() * scaldot(du2,v1) );
  for(int i=0;i<3;i++) dt[i]= dt[i]/lb->getJacobian();
}


static void compute_Deltat_tilde(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int &n, double Deltat[256][3][3]){
  SVector3 vect(0.,0.,0.),vect2(0.,0.,0.);
  vect = crossprod(lb->gett0(),lb->getphi0(0));
  vect2= crossprod(lb->gett0(),lb->getphi0(1));
  double invJ0 = 1./lb->getJacobian();
  for(int i=0;i<n;i++){
    Deltat[i][0][0] = invJ0 *(-lb->gett0(0)*vect(2)*Grads[i](1) + lb->gett0(0)*vect2(2))*Grads[i](0);
    Deltat[i][1][0] = invJ0 *((lb->getphi0(0,2)-lb->gett0(1)*vect(2))*Grads[i](1) - (lb->getphi0(1,2)-lb->gett0(1)*vect2(2))*Grads[i](0));
    Deltat[i][2][0] = invJ0 *((-lb->getphi0(0,1)-lb->gett0(2)*vect(2))*Grads[i](1) - (-lb->getphi0(1,1)-lb->gett0(2)*vect2(2))*Grads[i](0));

    Deltat[i][0][1] = invJ0 *((-lb->getphi0(0,2)-lb->gett0(0)*vect(2))*Grads[i](1) - (-lb->getphi0(1,2)-lb->gett0(0)*vect2(2))*Grads[i](0));
    Deltat[i][1][1] = invJ0 *(-lb->gett0(1)*vect(2)*Grads[i](1) + lb->gett0(1)*vect2(2)*Grads[i](0));
    Deltat[i][2][1] = invJ0 *((lb->getphi0(0,0)-lb->gett0(2)*vect(2))*Grads[i](1) - (lb->getphi0(1,0)-lb->gett0(2)*vect2(2))*Grads[i](0));

    Deltat[i][0][2] = invJ0 *((lb->getphi0(0,1)-lb->gett0(0)*vect(2))*Grads[i](1) - (lb->getphi0(1,1)-lb->gett0(0)*vect2(2))*Grads[i](0));
    Deltat[i][1][2] = invJ0 *((-lb->getphi0(0,0)-lb->gett0(1)*vect(2))*Grads[i](1) - (-lb->getphi0(1,0)-lb->gett0(1)*vect2(2))*Grads[i](0));
    Deltat[i][2][2] = invJ0 *( -lb->gett0(2)*vect(2)*Grads[i](1) + lb->gett0(2)*vect2(2)*Grads[i](0));
  }
}

static void matTvectprod(const double m[3][3], const SVector3 &v, double v2[3]){
  double temp[3];
  for(int i=0;i<3;i++){
    temp[0] = m[0][i]; temp[1] = m[1][i]; temp[2] = m[2][i];
    v2[i] = dot(temp,v);
  }
}

static void matvectprod(const double m[3][3], const SVector3 &v, double v2[3]){
  double temp[3];
  for(int i=0;i<3;i++){
    temp[0] = m[i][0]; temp[1] = m[i][1]; temp[2] = m[i][2];
    v2[i] = dot(temp,v);
  }
}

static inline double scaldot(const SVector3 &a, const SVector3 &b){
  double c=0.;
  for(int i=0;i<3;i++)
    c+= a[i]*b[i];
  return c;
}
static inline void diff(const SVector3 &a,const SVector3 &b,SVector3 &c){
  for(int i=0;i<3;i++)
    c[i]=a[i]-b[i];
}

double epsilongd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::GradType> &Grads, const std::vector<double> &disp){
  int nbFF=Grads.size();
  double eps = 0.,u_g,u_d;
  SVector3 u;
  for(int i=0;i<nbFF;i++){
    u(0)=disp[i]; u(1)=disp[i+nbFF]; u(2)=disp[i+2*nbFF];
    u_g=scaldot(u,lb->getphi0(gamma));
    u_d=scaldot(u,lb->getphi0(delta));
    eps+=(u_g*Grads[i](delta)+u_d*Grads[i](gamma));
  }
  eps/=2.;
  return eps;
}

double rhogd(const int gamma, const int delta, const LocalBasis *lb,const std::vector<TensorialTraits<double>::HessType> &hess, const std::vector<double> &disp){
  int nbFF=hess.size();
  double rho = 0.,u_t;
  SVector3 u;
  for(int i=0;i<nbFF;i++){
    u(0)=disp[i]; u(1)=disp[i+nbFF]; u(2)=disp[i+2*nbFF];
    u_t=scaldot(u,lb->gett0());
    rho-=(u_t*hess[i](gamma,delta));
  }
  return rho;
}

void stressReduction(const LinearElasticShellHookeTensor *H,const std::vector<TensorialTraits<double>::GradType> &Grads,
                     const LocalBasis *lb,const std::vector<double> &disp, reductionElement &n){
  n.setAll(0.);
  double eps_gd;
  for(int gamma=0;gamma<2;gamma++)
    for(int delta=0;delta<2;delta++)
    {
      eps_gd=epsilongd(gamma,delta,lb,Grads,disp);
        for(int alpha=0;alpha<2;alpha++)
          for(int beta=0;beta<2;beta++)
            n(alpha,beta) += H->get(alpha,beta,gamma,delta)*eps_gd;
    }
}

void momentReduction(const LinearElasticShellHookeTensor *H,const std::vector<TensorialTraits<double>::HessType> &hess,
                     const LocalBasis *lb,const std::vector<double> &disp, reductionElement &m){
  m.setAll(0.);
  double rho_gd;
  for(int gamma=0;gamma<2;gamma++)
    for(int delta=0;delta<2;delta++)
    {
      rho_gd=rhogd(gamma,delta,lb,hess,disp);
      for(int alpha=0;alpha<2;alpha++)
        for(int beta=0;beta<2;beta++)
          m(alpha,beta) += H->get(alpha,beta,gamma,delta)*rho_gd;
    }
}

void stressReductionHat(const reductionElement &n,const LocalBasis *lb, reductionElement &nhat){
  nhat.setAll(0.);
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++){
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++)
          nhat(alpha,beta)+=lb->getT(gamma,alpha)*n(gamma,delta)*lb->getT(delta,beta);
    }
}
// Should be somewhere else ??
void displacementjump(const std::vector<double> &Val_m,const int n_m,const std::vector<double> &Val_p,const int n_p,const std::vector<double> &disp,SVector3 &ujump){
  SVector3 up(0.,0.,0.),um(0.,0.,0.);
  for(int i=0;i<3;i++){
    for(int j=0;j<n_m;j++)
      um(i)+=Val_m[j]*disp[j+i*n_m];
    for(int j=0;j<n_p;j++)
      up(i)+=Val_p[j]*disp[j+i*n_p+3*n_m];
  }
  ujump = up-um;
}

void displacementjump(const std::vector<double> &Val_m,const int n_m,const std::vector<double> &Val_p,
                      const int n_p,const std::vector<double> &dispm,
                      const std::vector<double> &dispp,SVector3 &ujump){
  SVector3 up(0.,0.,0.),um(0.,0.,0.);
  for(int i=0;i<3;i++){
    for(int j=0;j<n_m;j++)
      um(i)+=Val_m[j]*dispm[j+i*n_m];
    for(int j=0;j<n_p;j++)
      up(i)+=Val_p[j]*dispp[j+i*n_p];
  }
  ujump = up-um;
}


void rotationjump(const LocalBasis *lbs,const std::vector<SVector3> &Gradm,const int nbFFm, const LocalBasis *lbm,
                  const std::vector<SVector3> &Gradp,const int nbFFp,const LocalBasis *lbp,
                  const std::vector<double> &dispm, const std::vector<double> &dispp, double rjump[3])
{
  for(int i=0;i<3;i++) rjump[i] =0.;
  // Compute Dt linear case
  /*double Deltat[256][3][3];
  double dt[3];
  SVector3 u(0.,0.,0.);
  compute_Deltat_tilde(lbm, Gradm, nbFFm, Deltat);
  for(int i=0;i<nbFFm;i++){
    u[0] = dispm[i]; u[1] = dispm[i+nbFFm]; u[2]=dispm[i+nbFFm+nbFFm];
    matvectprod(Deltat[i],u,dt);
    for(int j=0;j<3;j++)
      rjump[j]-= dt[j];
  }
  compute_Deltat_tilde(lbp, Gradp, nbFFp, Deltat);
  for(int i=0;i<nbFFp;i++){
    u[0] = dispp[i]; u[1] = dispp[i+nbFFp]; u[2] = dispp[i+nbFFp+nbFFp];
    matvectprod(Deltat[i],u,dt);
    for(int j=0;j<3;j++)
      rjump[j] += dt[j];
  }*/
  // general case
  SVector3 dtt;
  SVector3 dt;
  SVector3 dtm = lbm->gett0();
  SVector3 dtp = lbp->gett0();
  normalVariation(lbm,nbFFm,Gradm,dispm,dtt);
  for(int i=0;i<3;i++)
    dtm[i] += dtt[i];
  normalVariation(lbp,nbFFp,Gradp,dispp,dtt);
  for(int i=0;i<3;i++)
    dtp[i] += dtt[i];
  dt = dtp -dtm;
  // dt is a vector with normal variation
  //SVector3 t;
  //t=lbs->gett0() +dt;
  //rjump[0] = - atan(scaldot(dt,lbs->getphi0(0))/(abs(lbs->getphi0(0,1))*scaldot(dt,lbs->getphi0(1))));
  //rjump[1] = - atan(abs(lbs->getphi0(0,1))*scaldot(dt,lbs->getphi0(1))/scaldot(dt,lbs->getphi0(0))); // As phi0,0 is not normed
  for(int alpha=0;alpha<2;alpha++)
    rjump[alpha] = scaldot(dt,lbs->getphi0d(alpha));
}

void rotationjump(const LocalBasis *lbs,const std::vector<SVector3> &Gradm,const int nbFFm, const int nbdofm,const LocalBasis *lbm,
                  const std::vector<SVector3> &Gradp,const int nbFFp,const LocalBasis *lbp,
                  const std::vector<double> &disp,  double rjump[3])
{
/*  // Compute Dt
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) rjump[i][j] =0.;
  double Deltat[256][3][3];
  double dt[3];
  compute_Deltat_tilde(lbm, Gradm, nbFFm, Deltat);
  for(int i=0;i<nbFFm;i++)
    for(int alpha=0;alpha<2;alpha++){
      matTvectprod(Deltat[i],lbm->getphi0(alpha),dt);
      for(int j=0;j<3;j++)
        rjump[j][alpha] -= dt[j]*disp[i+j*nbFFm];
    }
  compute_Deltat_tilde(lbp, Gradp, nbFFp, Deltat);
  for(int i=0;i<nbFFp;i++)
    for(int alpha=0;alpha<2;alpha++){
      matTvectprod(Deltat[i],lbp->getphi0(alpha),dt);
      for(int j=0;j<3;j++)
        rjump[j][alpha] += dt[j]*disp[i+j*nbFFp+nbdofm];
    }*/
//  for(int i=0;i<3;i++) rjump[i] =0.;
/*  // Compute Dt linear case
  double Deltat[256][3][3];
  double dt[3];
  SVector3 u(0.,0.,0.);
  compute_Deltat_tilde(lbm, Gradm, nbFFm, Deltat);
  for(int i=0;i<nbFFm;i++){
    u[0] = disp[i]; u[1] = disp[i+nbFFm]; u[2] = disp[i+nbFFm+nbFFm];
    matvectprod(Deltat[i],u,dt);
    for(int j=0;j<3;j++)
      rjump[j] -= dt[j];
  }
  compute_Deltat_tilde(lbp, Gradp, nbFFp, Deltat);
  for(int i=0;i<nbFFp;i++){
    u[0] = disp[i+nbdofm]; u[1] = disp[i+nbFFp+nbdofm]; u[2] = disp[i+nbFFp+nbFFp+nbdofm];
    matvectprod(Deltat[i],u,dt);
    for(int j=0;j<3;j++)
      rjump[j] += dt[j];
  }*/
  // general case
  SVector3 dtt;
  SVector3 dtm = lbm->gett0();
  SVector3 dtp = lbp->gett0();
  SVector3 dt;
  normalVariation(lbm,nbFFm,Gradm,disp,dtt);
  for(int i=0;i<3;i++)
    dtm[i] += dtt[i];
  normalVariation(lbp,nbFFp,Gradp,disp,dtt,nbdofm);
  for(int i=0;i<3;i++)
    dtp[i] += dtt[i];
  dt = dtp-dtm;

  // dt is a vector with normal variation
  // now projection to have the angular component
  //rjump[0] = - atan(scaldot(dt,lbs->getphi0(0))/(abs(lbs->getphi0(0,1))*scaldot(dt,lbs->getphi0(1))));
  //rjump[1] = - atan(abs(lbs->getphi0(0,1))*scaldot(dt,lbs->getphi0(1))/scaldot(dt,lbs->getphi0(0))); // As phi0,0 is not normed
  for(int alpha=0;alpha<2;alpha++)
    rjump[alpha] = scaldot(dt,lbs->getphi0d(alpha));
}
