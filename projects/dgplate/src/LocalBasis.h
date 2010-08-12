//
// C++ Interface: terms
//
// Description: Define the localbasis of element for shells
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef LOCALBASIS_H_
#define LOCALBASIS_H_
#include "SVector3.h"
#include "MElement.h"
#include "functionSpace.h"
#include "MInterfaceElement.h"
class LocalBasis{
  protected :
    std::vector<SVector3> phi0;
    std::vector<SVector3> phi0d;
    SVector3 derivphi0[2][2];
    fullMatrix<double> T; // not used if interface element  and bulk term --> create separate class ??
    fullMatrix<double> t; // not used if interface element  and bulk term --> create separate class ??
    SVector3 t0;
    double detJ0;

  void setphiall(const double d){
    for(int i=0;i<2;i++)
      for(int j=0;j<3;j++){
        phi0[i](j)=d;
        phi0d[i](j)=d;
        for(int jj=0;jj<2;jj++)
          derivphi0[i][jj](j) = d;
      }
  }

  public :
    LocalBasis(){
      phi0.reserve(2);
      phi0d.reserve(2);
      T.resize(2,2);
      t.resize(2,2);
      SVector3 temp(0.,0.,0.);
      for(int i=0;i<2;i++){
        phi0.push_back(temp);
        phi0d.push_back(temp);
        for(int j=0;j<2;j++){
          derivphi0[i][j]=temp;
        }
      }
      T.setAll(0.);
      t.setAll(0.);
    }
    ~LocalBasis(){}
    // Copy Constructor
    LocalBasis& operator=(const LocalBasis &source){
      phi0 = source.phi0;
      phi0d = source.phi0d;
      T = source.T;
      t = source.t;
      t0= source.t0;
      detJ0 = source.detJ0;
      for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
          derivphi0[i][j] = source.derivphi0[i][j];
    }
    LocalBasis(const LocalBasis &source){
      phi0 = source.phi0;
      phi0d = source.phi0d;
      T = source.T;
      t = source.t;
      t0= source.t0;
      detJ0 = source.detJ0;
      for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
          derivphi0[i][j] = source.derivphi0[i][j];
    }

    void set(MElement *ele, const std::vector<TensorialTraits<double>::GradType> &Grads,
             const std::vector<TensorialTraits<double>::HessType> & Hess){
      const int nbFF = ele->getNumVertices();

      // Compute local basis vector
      this->setphiall(0.);
      for(int i=0;i<2;i++){
        for(int j=0;j<nbFF;j++){
          phi0[i](0) += Grads[j](i)*ele->getVertex(j)->x();
          phi0[i](1) += Grads[j](i)*ele->getVertex(j)->y();
          phi0[i](2) += Grads[j](i)*ele->getVertex(j)->z();

          for(int jj=0;jj<2;jj++){
            derivphi0[i][jj](0) += Hess[j](i,jj)*ele->getVertex(j)->x();
            derivphi0[i][jj](1) += Hess[j](i,jj)*ele->getVertex(j)->y();
            derivphi0[i][jj](2) += Hess[j](i,jj)*ele->getVertex(j)->z();
          }
        }
      }

      // Compute normal and Jacobian
      t0 = crossprod(phi0[0],phi0[1]);
      detJ0=t0.norm();
      t0.normalize();

      // Compute dual basis vector
      double cos;
      for(int i=0;i<2;i++){
        phi0d[i] = crossprod(phi0[1-i],t0);
        cos =  dot(phi0d[i],phi0[i]);
        phi0d[i] *= (1./cos);
      }
    }

    void set(MInterfaceElement *iele, const std::vector<TensorialTraits<double>::GradType> &Grads, const SVector3 &t0p, const SVector3 &t0m){
      const int nbFF = iele->getNumVertices();

      // Compute of normal
      t0=t0p+t0m;
      //assert(t0.norm()>1.e-8);
      t0.normalize();

      // Compute local basis vector
      this->setphiall(0.);
      for(int j=0;j<nbFF;j++){
          phi0[0](0) += Grads[j](0)*iele->getVertex(j)->x();
          phi0[0](1) += Grads[j](0)*iele->getVertex(j)->y();
          phi0[0](2) += Grads[j](0)*iele->getVertex(j)->z();
        }
        phi0[1] = crossprod(t0,phi0[0]);
        phi0[1].normalize();

      // Compute normal and Jacobian
      detJ0=crossprod(phi0[0],phi0[1]).norm();

      // Compute dual basis vector
      double cos;
      for(int i=0;i<2;i++){
        phi0d[i] = crossprod(phi0[1-i],t0);
        cos =  dot(phi0d[i],phi0[i]);
        phi0d[i] *= (1./cos);
      }
    }

    void set(MInterfaceElement *iele, const std::vector<TensorialTraits<double>::GradType> &Grads, const SVector3 &t0m){
      const int nbFF = iele->getNumVertices();

      // Compute of normal
      t0=t0m;
      //assert(t0.norm()>1.e-8);
      t0.normalize();

      // Compute local basis vector
      this->setphiall(0.);
      for(int j=0;j<nbFF;j++){
          phi0[0](0) += Grads[j](0)*iele->getVertex(j)->x();
          phi0[0](1) += Grads[j](0)*iele->getVertex(j)->y();
          phi0[0](2) += Grads[j](0)*iele->getVertex(j)->z();
        }
        phi0[1] = crossprod(t0,phi0[0]);
        phi0[1].normalize();

      // Compute normal and Jacobian
      detJ0=crossprod(phi0[0],phi0[1]).norm();

      // Compute dual basis vector
      double cos;
      for(int i=0;i<2;i++){
        phi0d[i] = crossprod(phi0[1-i],t0);
        cos =  dot(phi0d[i],phi0[i]);
        phi0d[i] *= (1./cos);
      }
    }


    // Push Forward Tensor (not used for interface element --> create a separate class ??
    void set_pushForward(const LocalBasis *lbs){
      // push forward tensor T and inverse push forward tensor t
      for(int i=0;i<2;i++)
        for(int j=0;j<2;j++){
          T(i,j)=dot(phi0[i],lbs->getphi0d(j));
          t(i,j)=dot(phi0d[j],lbs->getphi0(i));
        }
    }
    // get operation
    inline double getJacobian() const {return detJ0;}
    SVector3 gett0() const {return t0;}
    inline double gett0(const int i) const {return t0(i);}
    std::vector<SVector3> getphi0() const {return phi0;}
    std::vector<SVector3> getphi0d()const {return phi0d;}
    const SVector3& getphi0(const int i) const {if(i<2) return phi0[i]; else return t0;}
    const SVector3& getphi0d(const int i)const {if (i<2) return phi0d[i]; else return t0;}
    inline double getphi0(const int i,const int j) const {return phi0[i](j);}
    inline double getphi0d(const int i,const int j) const {return phi0d[i](j);}
    fullMatrix<double> getT() const {return T;}
    fullMatrix<double> gett() const {return t;}
    inline double getT(const int i, const int j) const {return T(i,j);}
    inline double gett(const int i, const int j) const {return t(i,j);}
    // get a orthonormal basis
    inline SVector3 getOrthonormalVector(const int i) const{
    if(i<2)
      return phi0[i]*(1/phi0[i].norm());
    else
      return t0*(1/t0.norm());
    }
    inline const SVector3& getderivphi0(const int i, const int j) const {return derivphi0[i][j];}
    inline const double getderivphi0(const int i, const int j, const int k) const {return derivphi0[i][j](k);}

    // Print
    void print() const {
    printf("Basis phi0\n");
    printf("1 : %f %f %f\n",phi0[0](0),phi0[0](1),phi0[0](2));
    printf("2 : %f %f %f\n",phi0[1](0),phi0[1](1),phi0[1](2));
    printf("Normal : %f %f %f \n",t0(0),t0(1),t0(2));
    printf("Dual Basis phi0d\n");
    printf("1 : %f %f %f\n",phi0d[0](0),phi0d[0](1),phi0d[0](2));
    printf("2 : %f %f %f\n",phi0d[1](0),phi0d[1](1),phi0d[1](2));
    printf("Derivative of phi0\n");
    printf("11 : %lf %lf %lf\n", derivphi0[0][0](0),derivphi0[0][0](1),derivphi0[0][0](2));
    printf("12 : %lf %lf %lf\n", derivphi0[0][1](0),derivphi0[0][1](1),derivphi0[0][1](2));
    printf("21 : %lf %lf %lf\n", derivphi0[1][0](0),derivphi0[1][0](1),derivphi0[1][0](2));
    printf("22 : %lf %lf %lf\n", derivphi0[1][1](0),derivphi0[1][1](1),derivphi0[1][1](2));
    T.print("PushForward tensor");
    t.print("Inverse PushForward tensor");
    }
};
# endif // localbasis
