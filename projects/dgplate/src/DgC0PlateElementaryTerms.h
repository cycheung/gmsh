//
// C++ Interface: terms
//
// Description: Functions used to compute a component of elementary matrix term
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef DGC0PLATEELEMENTARYTERMS_H_
#define DGC0PLATEELEMENTARYTERMS_H_
#include "SVector3.h"
#include "functionSpace.h"
#include "LinearElasticShellHookeTensor.h"
inline void diaprod(const double a[3], const double b[3], double m[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      m[i][j]=a[i]*b[j];
}

inline void diaprodAdd(const double a[3], const double b[3], double m[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      m[i][j]+=a[i]*b[j];
}

inline void diaprod(const double a[3], const SVector3 &b, double m[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      m[i][j]=a[i]*b[j];
}

inline void diaprod(const SVector3 &a, const double b[3], double m[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      m[i][j]=a[i]*b[j];
}

inline void diaprod(const SVector3 &a, const SVector3 &b, double m[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      m[i][j]=a[i]*b[j];
}

inline double scaldot(const double a[3],const SVector3 &b){
  double c=0.;
  for(int i=0;i<3;i++)
    c+= a[i]*b[i];
  return c;
}

static inline double scaldot(const SVector3 &a,const SVector3 &b){
  double c=0.;
  for(int i=0;i<3;i++)
    c+= a[i]*b[i];
  return c;
}

static inline void dot(const double a[3], const double b[3], double c[3]){
  for(int i=0;i<3;i++)
    c[i] = a[i]*b[i];
}


static inline void dot(const double a[3], const SVector3 &b, double c[3]){
  for(int i=0;i<3;i++)
    c[i] = a[i]*b(i);
}

inline void matTvectprod(const double m[3][3], const SVector3 &v, double v2[3]){
  double temp[3];
  for(int i=0;i<3;i++){
    temp[0] = m[0][i]; temp[1] = m[1][i]; temp[2] = m[2][i];
    v2[i] = dot(temp,v);
  }
}
inline void matTvectprod(const double m[3][3], const SVector3 &v, SVector3 &v2){
  v2[0] = 0.; v2[1]=0.; v2[2]=0.;
  double temp[3];
  for(int i=0;i<3;i++){
    temp[0] = m[0][i]; temp[1] = m[1][i]; temp[2] = m[2][i];
    v2[i] = dot(temp,v);
  }
}


inline void matTvectprod(const double m[3][3], const double v[3],  double v2[3]){
  double temp[3];
  for(int i=0;i<3;i++){
    temp[0] = m[0][i]; temp[1] = m[1][i]; temp[2] = m[2][i];
    v2[i] = dot(temp,v);
  }
}
static inline void matvectprod(const double m[3][3], const SVector3 &v1,double v2[3]){
  for(int i=0;i<3;i++)
    v2[i]=scaldot(m[i],v1);
}

// Compute the component (j,k) of the elementary stiffness matrix
inline void BulkC0PlateDGStiffnessMembraneTerms(const double Bj[3][2][2],const double Bk[3][2][2], const LinearElasticShellHookeTensor *H, double me[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      me[i][j]=0.; // make a methode in object me ??
  double mtemp[3][3];
  double v1[3], v2[3];
  for(int a=0;a<2;a++)
    for(int b=0;b<2;b++){
      v2[0] = Bk[0][a][b]; v2[1] = Bk[1][a][b]; v2[2]= Bk[2][a][b]; // make a method
      for(int c=0;c<2;c++)
        for(int d=0;d<2;d++){
          v1[0] = Bj[0][c][d]; v1[1]=Bj[1][c][d]; v1[2]=Bj[2][c][d]; //make a method
          diaprod(v1,v2,mtemp);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              me[jj][kk]+=(H->get(a,b,c,d)*mtemp[jj][kk]);
        }
    }
}

inline void BulkC0PlateDGStiffnessBendingTerms(const double Bmj[3][2][2], const double Bmk[3][2][2], const LinearElasticShellHookeTensor *H, double me[3][3]){
  double mtemp[3][3];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double B1[3],B2[3];
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++){
      B1[0]= Bmk[0][alpha][beta]; B1[1]= Bmk[1][alpha][beta]; B1[2]= Bmk[2][alpha][beta];
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          B2[0]= Bmj[0][gamma][delta]; B2[1]= Bmj[1][gamma][delta]; B2[2]= Bmj[2][gamma][delta];
          diaprod(B2,B1,mtemp);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              me[jj][kk]+=H->get(alpha,beta,gamma,delta)*mtemp[jj][kk];
        }
    }
}

inline void consC0PlateStiffnessTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double dt[3][3], const LocalBasis *lb, double me[3][3]){
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double v1[3], v2[3];
  double temp=0.;
  double temp2=0.;
  for(int alpha=0;alpha<2;alpha++){
    matTvectprod(dt,lb->getphi0(alpha),v2);
    for(int beta=0;beta<2;beta++){
      temp = 0.5*(-scaldot(lb->getphi0(1),lb->getphi0(beta)));
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          temp2 = temp *Hhat->get(alpha,beta,gamma,delta);
          v1[0] = temp2*Bhat[0][gamma][delta]; v1[1] = temp2*Bhat[1][gamma][delta]; v1[2] = temp2*Bhat[2][gamma][delta];
          diaprodAdd(v2,v1,me);
        }
    }
  }
}

inline void compC0PlateStiffnessTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double dt[3][3], const LocalBasis *lb, double me[3][3]){
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double v1[3], v2[3];
  double mtemp[3][3];
  double temp=0.;
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++)
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          v1[0] = Bhat[0][gamma][delta]; v1[1] = Bhat[1][gamma][delta]; v1[2] = Bhat[2][gamma][delta];
          matTvectprod(dt,lb->getphi0(alpha),v2);
          diaprod(v1,v2,mtemp);
          temp = 0.5*Hhat->get(alpha,beta,gamma,delta)*(-scaldot(lb->getphi0(1),lb->getphi0(beta)));
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
              me[i][j] += (mtemp[i][j]*temp);
        }
}

inline void stabilityC0PlateStiffnessTerms(LinearElasticShellHookeTensor *Hhat, const double dta[3][3], const double dtb[3][3], const LocalBasis *lb, double me[3][3]){
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double mtemp[3][3];
  double v1[3],v2[3];
  double temp=0.;
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++)
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          matTvectprod(dta,lb->getphi0(gamma),v1);
          matTvectprod(dtb,lb->getphi0(alpha),v2);
          diaprod(v1,v2,mtemp);
          temp = Hhat->get(alpha,beta,gamma,delta)*(-scaldot(lb->getphi0(1),lb->getphi0(delta)))*(-scaldot(lb->getphi0(1),lb->getphi0(beta)));
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
              me[i][j]+=(mtemp[i][j]*temp);
        }
}

inline void compC0PlateStiffnessMembraneTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double Na, const LocalBasis *lb, double me[3][3]){
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double v1[3];// v2[3];
  double mtemp[3][3];
  double temp=0.;
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++)
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          v1[0] = Bhat[0][gamma][delta]; v1[1] = Bhat[1][gamma][delta]; v1[2] = Bhat[2][gamma][delta];
          diaprod(v1,lb->getphi0(beta),mtemp);
          temp = 0.5*Na*Hhat->get(alpha,beta,gamma,delta)*(-scaldot(lb->getphi0(1),lb->getphi0(alpha)));
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
              me[i][j] += (mtemp[i][j]*temp);
        }
}

inline void compC0PlateForceMembraneTerms(const int beta,const int gamma,const int delta,const LocalBasis *lb, const LocalBasis *lbs, const TensorialTraits<double>::GradType &Gradj,const SVector3 &ujump,double me_comp[3]){
  for(int i=0;i<3;i++) me_comp[i]=0.;
  SVector3 v4(0.,0.,0.);
  double v1;
  double mtemp[3][3];
  for(int i=0;i<3;i++){
    for(int alpha=0;alpha<2;alpha++)
      for(int mu=0;mu<2;mu++){
        v4(i)+= lb->gett(gamma,alpha)*(0.5*(Gradj(alpha)*lb->getphi0(mu,i)+Gradj(mu)*lb->getphi0(alpha,i)))*lb->gett(delta,mu);
      }
  }
  diaprod(v4,lbs->getphi0(beta),mtemp); // normallement l'inverse ??
  matvectprod(mtemp,ujump,me_comp);
}

inline void consC0PlateStiffnessMembraneTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double Na, const LocalBasis *lb, double me[3][3]){
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double v1[3];
  double mtemp[3][3];
  double temp=0.;
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++)
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          temp = 0.5*Na*Hhat->get(alpha,beta,gamma,delta)*(-scaldot(lb->getphi0(1),lb->getphi0(alpha)));
          v1[0] = temp*Bhat[0][gamma][delta]; v1[1] = temp*Bhat[1][gamma][delta]; v1[2] = temp*Bhat[2][gamma][delta];
          diaprod(lb->getphi0(beta),v1,mtemp);
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
              me[i][j] += (mtemp[i][j]);
        }
}

inline void stabilityC0PlateStiffnessMembraneTerms(LinearElasticShellHookeTensor *Hhat, const double Na, const double Nb, const LocalBasis *lb, double me[3][3]){
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) me[i][j]=0.;
  double mtemp[3][3];
  double temp=0.;
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++)
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          diaprod(lb->getphi0(gamma),lb->getphi0(alpha),mtemp);
          temp = Na*Nb*Hhat->get(alpha,beta,gamma,delta)*(-scaldot(lb->getphi0(1),lb->getphi0(delta)))*(-scaldot(lb->getphi0(1),lb->getphi0(beta)));
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
              me[i][j]+=(mtemp[i][j]*temp);
        }
}

inline void stabilityC0PlateStiffnessShearingTerms(const double Bj[3], const double Bk[3],double me[3][3]){
  diaprod(Bj,Bk,me);
}

inline void stabilityC0PlateForceMembraneTerms(LinearElasticShellHookeTensor *Hhat, const std::vector<SVector3> &Vals_m,const std::vector<SVector3> &Vals_p, const double Nb, const LocalBasis *lb, const std::vector<double> &disp,double me[3]){
  for(int i=0;i<3;i++) me[i]=0.;
  double mtemp[3][3];
  double temp;
  int n_m=Vals_m.size();
  int n_p=Vals_p.size();
  for(int alpha=0;alpha<2;alpha++)
    for(int beta=0;beta<2;beta++)
      for(int gamma=0;gamma<2;gamma++)
        for(int delta=0;delta<2;delta++){
          diaprod(lb->getphi0(gamma),lb->getphi0(alpha),mtemp);
          for(int j=0;j<n_m;j++){
            temp = - Vals_m[j](0)*Nb*Hhat->get(alpha,beta,gamma,delta)*(-lb->getphi0(1,delta))*(-lb->getphi0(1,beta));
            for(int i=0;i<3;i++)
              me[i]+=temp*(mtemp[i][0]*disp[j]+mtemp[i][1]*disp[j+n_m]+mtemp[i][2]*disp[j+2*n_m]);
          }
          for(int j=0;j<n_p;j++){
            temp = Vals_p[j](0)*Nb*Hhat->get(alpha,beta,gamma,delta)*(-lb->getphi0(1,delta))*(-lb->getphi0(1,beta));
            for(int i=0;i<3;i++)
              me[i]+=temp*(mtemp[i][0]*disp[j+3*n_m]+mtemp[i][1]*disp[j+n_p+3*n_m]+mtemp[i][2]*disp[j+2*n_p+3*n_m]);
          }
        }
}

inline void stabC0PlateForceMembraneTerms(const int beta,const int delta,const LocalBasis *lb,const SVector3 &ujump,double me[3]){
  for(int i=0;i<3;i++) me[i]=0.;
  double mtemp[3][3];
  diaprod(lb->getphi0(beta),lb->getphi0(delta),mtemp);
  matvectprod(mtemp,ujump,me);
}

inline void stabilityC0PlateForceShearingTerms(const double Bj[3], const double B_m[256][3],const double B_p[256][3],const int n_m,const int n_p,const std::vector<double> &disp,double me[3]){
  me[0]=0.;me[1]=0.;me[2]=0.;
  double mtemp[3][3];
  for(int j=0;j<n_m;j++){
    diaprod(Bj,B_m[j],mtemp);
    for(int i=0;i<3;i++)
      me[i]+=-(mtemp[i][0]*disp[j]+mtemp[i][1]*disp[j+n_m]+mtemp[i][2]*disp[j+2*n_m]);
  }
  for(int j=0;j<n_p;j++){
    diaprod(Bj,B_p[j],mtemp);
    for(int i=0;i<3;i++)
      me[i]+=(mtemp[i][0]*disp[j+3*n_m]+mtemp[i][1]*disp[j+n_p+3*n_m]+mtemp[i][2]*disp[j+2*n_p+3*n_m]);
  }
}

inline void compC0PlateForceTerms(const LinearElasticShellHookeTensor *Hhat, const int n_m, const int n_p, const double Bhat[3][2][2], const double Dt_m[256][3][3], const double Dt_p[256][3][3],const LocalBasis *lb, const std::vector<double> &disp, double me[3]){
  double mtemp[3][3];
  double v1[3],v2[3];
  double temp=0.;
  for(int i=0;i<3;i++) me[i]=0.;
  for(int a=0;a<2;a++)
    for(int b=0;b<2;b++)
      for(int c=0;c<2;c++)
        for(int d=0;d<2;d++){
          v1[0]=Bhat[0][c][d] ; v1[1]=Bhat[1][c][d]; v1[2]=Bhat[2][c][d];
          temp = 0.5*Hhat->get(a,b,c,d)*(-scaldot(lb->getphi0(1),lb->getphi0(b)));
          for(int j=0;j<n_m;j++){
            matTvectprod(Dt_m[j],lb->getphi0(a),v2); // put in a loop but in this time a vector is need to store the result of matTvectprod
            diaprod(v1,v2,mtemp);
            for(int k=0;k<3;k++)
              me[k] += -temp*(mtemp[k][0]*disp[j]+mtemp[k][1]*disp[j+n_m]+mtemp[k][2]*disp[j+2*n_m]);
          }
          for(int j=0;j<n_p;j++){
            matTvectprod(Dt_p[j],lb->getphi0(a),v2); //idem
            diaprod(v1,v2,mtemp);
            for(int k=0;k<3;k++)
              me[k] += temp*(mtemp[k][0]*disp[j+3*n_m]+mtemp[k][1]*disp[j+n_p+3*n_m]+mtemp[k][2]*disp[j+2*n_p+3*n_m]);
          }
        }
}

inline void stabilityC0PlateForceTerms(const int n_p,const int n_m, const LinearElasticShellHookeTensor *Hhat,const  double Dt[3][3],const double Dt_m[256][3][3],const double Dt_p[256][3][3], const LocalBasis *lb, const std::vector<double> &disp, double me[3]){
  for(int i=0;i<3;i++) me[i]=0.;
  double v1[3],v2[3];
  double mtemp[3][3];
  double temp=0.;
  for(int a=0;a<2;a++)
    for(int b=0;b<2;b++)
      for(int c=0;c<2;c++){
        matTvectprod(Dt,lb->getphi0(c),v1);
        for(int d=0;d<2;d++){
          temp = Hhat->get(a,b,c,d)*(-scaldot(lb->getphi0(1),lb->getphi0(d)))*(-scaldot(lb->getphi0(1),lb->getphi0(b)));
          for(int j=0;j<n_m;j++){
            matTvectprod(Dt_m[j],lb->getphi0(a),v2); // can be put on a loop but it is necessary to create a vector to keep the result of matTvectprod
            diaprod(v1,v2,mtemp);
            for(int jj=0;jj<3;jj++)
              me[jj] += -temp*(mtemp[jj][0]*disp[j]+mtemp[jj][1]*disp[j+n_m]+mtemp[jj][2]*disp[j+2*n_m]);
          }
          for(int j=0;j<n_p;j++){
            matTvectprod(Dt_p[j],lb->getphi0(a),v2); // idem
            diaprod(v1,v2,mtemp);
            for(int jj=0;jj<3;jj++)
              me[jj] += temp*(mtemp[jj][0]*disp[j+3*n_m]+mtemp[jj][1]*disp[j+n_p+3*n_m]+mtemp[jj][2]*disp[j+2*n_p+3*n_m]);
          }
       }
      }
}

inline void compC0PlateForceTermsBound(const LinearElasticShellHookeTensor *Hhat, const int n, const double Bhat[3][2][2],
                                       const double Dt_m[256][3][3], const LocalBasis *lb, const std::vector<double> &disp,
                                       double me[3]){
  double mtemp[3][3];
  double v1[3],v2[3];
  double temp=0.;
  for(int i=0;i<3;i++) me[i]=0.;
  for(int a=0;a<2;a++)
    for(int b=0;b<2;b++)
      for(int c=0;c<2;c++)
        for(int d=0;d<2;d++){
          v1[0]=Bhat[0][c][d] ; v1[1]=Bhat[1][c][d]; v1[2]=Bhat[2][c][d];
          temp = 0.5*Hhat->get(a,b,c,d)*(-scaldot(lb->getphi0(1),lb->getphi0(b)));
          for(int j=0;j<n;j++){
            matTvectprod(Dt_m[j],lb->getphi0(a),v2); // put in a loop but in this time a vector is need to store the result of matTvectprod
            diaprod(v1,v2,mtemp);
            for(int k=0;k<3;k++)
              me[k] += -temp*(mtemp[k][0]*disp[j]+mtemp[k][1]*disp[j+n]+mtemp[k][2]*disp[j+2*n]);
          }
        }
}

inline void stabilityC0PlateForceTermsBound(const int n_m,const LinearElasticShellHookeTensor *Hhat,const double Dt[3][3],
                                            const double Dt_m[256][3][3], const LocalBasis *lb, const std::vector<double> &disp, double me[3]){
  for(int i=0;i<3;i++) me[i]=0.;
  double v1[3],v2[3];
  double mtemp[3][3];
  double temp=0.;
  for(int a=0;a<2;a++)
    for(int b=0;b<2;b++)
      for(int c=0;c<2;c++){
        matTvectprod(Dt,lb->getphi0(c),v1);
        for(int d=0;d<2;d++){
          temp = Hhat->get(a,b,c,d)*(-scaldot(lb->getphi0(1),lb->getphi0(d)))*(-scaldot(lb->getphi0(1),lb->getphi0(b)));
          for(int j=0;j<n_m;j++){
            matTvectprod(Dt_m[j],lb->getphi0(a),v2); // can be put on a loop but it is necessary to create a vector to keep the result of matTvectprod
            diaprod(v1,v2,mtemp);
            for(int jj=0;jj<3;jj++)
              me[jj] += -temp*(mtemp[jj][0]*disp[j]+mtemp[jj][1]*disp[j+n_m]+mtemp[jj][2]*disp[j+2*n_m]);
          }
        }
      }
}

// Compute value needed in get
void compute_Bn(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int n, double B[][3][2][2]){
  for(int i=0;i<n;i++){
    for(int a=0;a<2;a++)
      for(int b=0;b<2;b++){
        B[i][0][a][b] = 0.5*(lb->getphi0(a,0)*Grads[i](b)+lb->getphi0(b,0)*Grads[i](a) );
        B[i][1][a][b] = 0.5*(lb->getphi0(a,1)*Grads[i](b)+lb->getphi0(b,1)*Grads[i](a) );
        B[i][2][a][b] = 0.5*(lb->getphi0(a,2)*Grads[i](b)+lb->getphi0(b,2)*Grads[i](a) );
      }
  }

}

void  Compute_Bm(const LocalBasis *lb,const std::vector<TensorialTraits<double>::GradType> &Grads,
                 const std::vector<TensorialTraits<double>::HessType> &Hess, const int &n, double B[256][3][2][2]){
  SVector3 t1(0.,0.,0.), t2(0.,0.,0.), p1(0.,0.,0.), p2(0.,0.,0.);
  SVector3 t1t(0.,0.,0.),t2t(0.,0.,0.),p1t(0.,0.,0.),p2t(0.,0.,0.);
  double scalt,scalt1, scalt2, scalp1, scalp2;
  double invJ = 1./lb->getJacobian();
  for(int alpha=0;alpha<2;alpha++){
    for(int beta=0;beta<2;beta++){
      scalt = invJ * scaldot(lb->getderivphi0(alpha,beta),lb->gett0());
      t1 = crossprod(lb->getphi0(0),lb->gett0());
      t2 = crossprod(lb->getphi0(1),lb->gett0());
      p1 = crossprod(lb->getderivphi0(alpha,beta),lb->getphi0(0));
      p2 = crossprod(lb->getderivphi0(alpha,beta),lb->getphi0(1));
      for(int i=0;i<n;i++){
        scalt1 = scalt * Grads[i](0);
        scalt2 = scalt * Grads[i](1);
        scalp1 = invJ*Grads[i](0);
        scalp2 = invJ*Grads[i](1);
        t1t = scalt2 * t1;
        t2t = scalt1 * t2;
        p1t = scalp2 * p1;
        p2t = scalp1 * p2;
        t1t+=p1t;
        t2t+=p2t;
        t2t-=t1t;
        for(int ii=0;ii<3;ii++){
          B[i][ii][alpha][beta] = t2t(ii)-lb->gett0(ii)*Hess[i](alpha,beta);
        }
      }
    }
  }
}

// Compute value needed in getinter
void compute_Bs(const LocalBasis *lb, const std::vector<TensorialTraits<double>::ValType> &Vals, const int n, double B[][3]){
  for(int i=0;i<n;i++){
    for(int a=0;a<3;a++)
        B[i][a] = lb->gett0(a)*Vals[i];
      }
}

void compute_Bhat(const LocalBasis *lb, const int n, const double B[256][3][2][2], double Bhat[256][3][2][2]){
  for(int i=0;i<n;i++){
    for(int ii=0;ii<3;ii++){
      Bhat[i][ii][0][0] =0.; Bhat[i][ii][0][1]=0.; Bhat[i][ii][1][0]=0.; Bhat[i][ii][1][1]=0.;
      for(int j=0;j<2;j++)
        for(int k=0;k<2;k++){
          Bhat[i][ii][0][0] += lb->gett(0,j)*B[i][ii][j][k]*lb->gett(0,k);
          Bhat[i][ii][0][1] += lb->gett(0,j)*B[i][ii][j][k]*lb->gett(1,k);
          Bhat[i][ii][1][0] += lb->gett(1,j)*B[i][ii][j][k]*lb->gett(0,k);
          Bhat[i][ii][1][1] += lb->gett(1,j)*B[i][ii][j][k]*lb->gett(1,k);
        }
    }
  }
}

// Delta tilde
void compute_Deltat_tilde(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int &n, double Deltat[256][3][3]){
  SVector3 vect(0.,0.,0.),vect2(0.,0.,0.);
  vect = crossprod(lb->gett0(),lb->getphi0(0));
  vect2= crossprod(lb->gett0(),lb->getphi0(1));
  double invJ0 = 1./lb->getJacobian();
  STensor3 mat1;
  STensor3 mat2;
  tensprod(lb->gett0(),vect,mat1);
  tensprod(lb->gett0(),vect2,mat2);

  double scr1[3][3];
  double scr2[3][3];
  scr1[0][0] =      0.;           scr1[0][1] = - lb->getphi0(0,2);  scr1[0][2] = lb->getphi0(0,1);
  scr1[1][0] = lb->getphi0(0,2);  scr1[1][1] = 0.  ;                scr1[1][2] = -lb->getphi0(0,0);
  scr1[2][0] = -lb->getphi0(0,1); scr1[2][1] = lb->getphi0(0,0);    scr1[2][2] = 0.;

  scr2[0][0] =      0.;           scr2[0][1] = - lb->getphi0(1,2);  scr2[0][2] = lb->getphi0(1,1);
  scr2[1][0] = lb->getphi0(1,2);  scr2[1][1] = 0.  ;                scr2[1][2] = -lb->getphi0(1,0);
  scr2[2][0] = -lb->getphi0(1,1); scr2[2][1] = lb->getphi0(1,0);    scr2[2][2] = 0.;
  for(int i=0;i<n;i++){
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
        Deltat[i][j][k] = invJ0*(Grads[i](1)*(scr1[j][k]-mat1(j,k))-Grads[i](0)*(scr2[j][k]-mat2(j,k)));
  }
}
void compute_Deltat_tildeBound(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int &n, double Deltat[256][3][3], const LocalBasis *lbs){
  SVector3 vect(0.,0.,0.),vect2(0.,0.,0.),temp2(0.,0.,0.),temp3(0.,0.,0.),phi1norm(0.,0.,0.);
  vect = crossprod(lb->gett0(),lb->getphi0(0));
  vect2= crossprod(lb->gett0(),lb->getphi0(1));
  phi1norm=lbs->getphi0(0);
  phi1norm.normalize();
  double temp[3][3];
  double invJ0 = 1./lb->getJacobian();
  STensor3 mat1;
  STensor3 mat2;
  tensprod(lb->gett0(),vect,mat1);
  tensprod(lb->gett0(),vect2,mat2);

  double scr1[3][3];
  double scr2[3][3];
  scr1[0][0] =      0.;           scr1[0][1] = - lb->getphi0(0,2);  scr1[0][2] = lb->getphi0(0,1);
  scr1[1][0] = lb->getphi0(0,2);  scr1[1][1] = 0.  ;                scr1[1][2] = -lb->getphi0(0,0);
  scr1[2][0] = -lb->getphi0(0,1); scr1[2][1] = lb->getphi0(0,0);    scr1[2][2] = 0.;

  scr2[0][0] =      0.;           scr2[0][1] = - lb->getphi0(1,2);  scr2[0][2] = lb->getphi0(1,1);
  scr2[1][0] = lb->getphi0(1,2);  scr2[1][1] = 0.  ;                scr2[1][2] = -lb->getphi0(1,0);
  scr2[2][0] = -lb->getphi0(1,1); scr2[2][1] = lb->getphi0(1,0);    scr2[2][2] = 0.;
  for(int i=0;i<n;i++){
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
        temp[j][k] = invJ0*(Grads[i](1)*(scr1[j][k]-mat1(j,k))-Grads[i](0)*(scr2[j][k]-mat2(j,k)));
    // Projection of normal
    for(int j=0;j<3;j++){
      temp3(0)=temp[0][j]; temp3(1)=temp[1][j]; temp3(2)=temp[2][j];
      temp2=dot(temp3,phi1norm)*phi1norm;
      temp3-=temp2;
      Deltat[i][0][j]=temp3(0);Deltat[i][1][j]=temp3(1); Deltat[i][2][j]=temp3(2);
    }
  }
};

void computeFintFrac(const MInterfaceElement *iele, const std::vector<double> &disp, const LocalBasis *lbm,
                             const LocalBasis *lbp, const LocalBasis *lbs, const double weight, const int numgauss, const int npts,
                             const SolElementType::eltype elemtype, const IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf,
                             const int nbFF_m, const std::vector<TensorialTraits<double>::ValType> &Vals_m,
                             const std::vector<TensorialTraits<double>::GradType> &Grads_m, const int nbdof_m,
                             const int nbFF_p, const std::vector<TensorialTraits<double>::ValType> &Vals_p,
                             const std::vector<TensorialTraits<double>::GradType> &Grads_p, fullVector<double> &m){

  double Deltat_m[256][3][3];
  double Deltat_p[256][3][3];
  reductionElement nhatmean;
  reductionElement mhatmean;
  ipf->getReductionFracture(iele,numgauss,npts,elemtype,disp,nbFF_m,nbdof_m,Vals_m,Grads_m,nbFF_p,Vals_p,Grads_p,nhatmean,mhatmean);
  double wJ = weight*lbs->getJacobian();
  // Compute of Deltat_tilde
  compute_Deltat_tilde(lbp,Grads_p,nbFF_p,Deltat_p);
  compute_Deltat_tilde(lbm,Grads_m,nbFF_m,Deltat_m);

  // Assemblage
  double dt[3];
  for(int alpha=0;alpha<2;alpha++){
    SVector3 malpha = mhatmean(alpha,0)*lbs->getphi0(0) + mhatmean(alpha,1)*lbs->getphi0(1);
    SVector3 nalpha = nhatmean(alpha,0) * lbs->getphi0(0) + nhatmean(alpha,1)*lbs->getphi0(1);
    for(int j=0;j<nbFF_m;j++){
      matTvectprod(Deltat_m[j],malpha,dt);
      for(int k=0;k<3;k++){
        m(j+k*nbFF_m)-= wJ*(dt[k]+nalpha(k)*Vals_m[j])*(-scaldot(lbs->getphi0(1),lbs->getphi0(alpha)));
      }
    }
    for(int j=0;j<nbFF_p;j++){
      matTvectprod(Deltat_p[j],malpha,dt);
      for(int k=0;k<3;k++)
        m(j+k*nbFF_p+nbdof_m)+= wJ*(dt[k]+nalpha(k)*Vals_p[j])*(-scaldot(lbs->getphi0(1),lbs->getphi0(alpha)));
    }
  }
};
#endif //DGC0PLATEELEMENTARYTERMS_H_

