//
// C++ Interface: terms
//
// Description: Elementary matrix terms for C0 Dg Plate
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef DGC0PLATETERMS_H
#define DGC0PLATETERMS_H
#include "SVector3.h"
#include <string>
#include "LocalBasis.h"
#include "LinearElasticShellHookeTensor.h"
#include "reduction.h"
#include "DgC0PlateFunctionSpace.h"
#include "simpleFunction.h"
#include "displacementField.h"
#include "IPField.h"
#include "terms.h"

template<class T1,class T2> class DgC0BilinearTerm : public BilinearTermBase
{
 protected :
  DgC0FunctionSpace<T1>& space1;
  DgC0FunctionSpace<T2>& space2;
 public :
  DgC0BilinearTerm(DgC0FunctionSpace<T1>& space1_,DgC0FunctionSpace<T1>& space2_) : space1(space1_),space2(space2_) {}
  virtual ~DgC0BilinearTerm() {}
};


template<class T1> class DgC0LinearTerm : public LinearTermBase
{
 protected :
  DgC0FunctionSpace<T1>& space1;
  bool inverseSign;
 public :
  DgC0LinearTerm(DgC0FunctionSpace<T1>& space1_) : space1(space1_), inverseSign(false) {}
  virtual ~DgC0LinearTerm() {}
  void invSign(){inverseSign ? inverseSign = false : inverseSign = true;}
};

template<class T1> class DgC0LoadTerm : public DgC0LinearTerm<T1>
{
  simpleFunction<typename TensorialTraits<T1>::ValType> &Load;
 public :
  DgC0LoadTerm(DgC0FunctionSpace<T1>& space1_,simpleFunction<typename TensorialTraits<T1>::ValType> &Load_) :DgC0LinearTerm<T1>(space1_),Load(Load_) {}
  virtual ~DgC0LoadTerm() {}

  virtual void get(MElement *ele,int npts,IntPt *GP,fullVector<double> &m)
  {
    int nbdof=DgC0LinearTerm<T1>::space1.getNumKeys(ele);
    int nbFF=nbdof/3;
    double jac[3][3];
    m.resize(nbdof);
    m.scale(0.);
    for (int i = 0; i < npts; i++)
    {
      const double u = GP[i].pt[0];const double v = GP[i].pt[1];const double w = GP[i].pt[2];
      const double weight = GP[i].weight;const double detJ = ele->getJacobian(u, v, w, jac);
      std::vector<double> Vals;
      DgC0LinearTerm<T1>::space1.f(ele,u, v, w, Vals);
      SPoint3 p;
      ele->pnt(u, v, w, p);
      typename TensorialTraits<T1>::ValType load=Load(p.x(),p.y(),p.z());
      for (int j = 0; j < nbFF ; ++j)
        for(int k=0;k<3;k++)
        m(j+k*nbFF)+=(Vals[j]*load(k))*weight*detJ;

    }
  }
};

template<class T1> class DgC0PressureLoadTerm : public DgC0LinearTerm<T1>
{
  simpleFunction<typename TensorialTraits<T1>::ValType> &Load;
 public :
  DgC0PressureLoadTerm(DgC0FunctionSpace<T1>& space1_,simpleFunction<typename TensorialTraits<T1>::ValType> &Load_) :DgC0LinearTerm<T1>(space1_),Load(Load_) {}
  virtual ~DgC0PressureLoadTerm() {}

  virtual void get(MElement *ele,int npts,IntPt *GP,fullVector<double> &m)
  {
    int nbdof=DgC0LinearTerm<T1>::space1.getNumKeys(ele);
    int nbFF=nbdof/3;
    double jac[3][3];
    m.resize(nbdof);
    m.scale(0.);
    SVector3 press = Load(1.,1.,1.); // Use a SVector3 because it is defined like this in NeumannBC change this ?
    for (int i = 0; i < npts; i++)
    {
      const double u = GP[i].pt[0];const double v = GP[i].pt[1];const double w = GP[i].pt[2];
      const double weight = GP[i].weight;const double detJ = ele->getJacobian(u, v, w, jac);
      std::vector<double> Vals;
      std::vector<SVector3> Grads;
      std::vector<STensor3> Hess;
      DgC0LinearTerm<T1>::space1.fuvw(ele,u, v, w, Vals);
      for (int j = 0; j < nbFF ; ++j){
        LocalBasis lb;
        double xyz[3];
        double uvw[3];
        xyz[0] = ele->getVertex(j)->x();
        xyz[1] = ele->getVertex(j)->y();
        xyz[2] = ele->getVertex(j)->z();
        ele->xyz2uvw(xyz,uvw);
        DgC0LinearTerm<T1>::space1.gradfuvw(ele,uvw[0],uvw[1],0.,Grads);
        DgC0LinearTerm<T1>::space1.hessfuvw(ele,uvw[0],uvw[1],0.,Hess);
        lb.set(ele,Grads,Hess);
        for(int k=0;k<3;k++)
          m(j+k*nbFF)-=Vals[j]*press(0)*lb.gett0(k)*weight*detJ; //
        Grads.clear(); Hess.clear();
      }
      Vals.clear();
    }
  }
};

class IsotropicElasticForceBulkTermC0Plate : public DgC0LinearTerm<SVector3>
{
 protected :
  double E,nu,h,Cm,Cn;
  bool fullDg;
  displacementField *ufield;
  IPField *ipf;
  SolElementType::eltype _elemtype;

 public :
  IsotropicElasticForceBulkTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,materialLaw* mlaw,double h_, bool FullDG,
                                       displacementField *uf,IPField *ip,
                                       SolElementType::eltype et) : DgC0LinearTerm<SVector3>(space1_),h(h_),
                                                                                   fullDg(FullDG),ufield(uf), ipf(ip),
                                                                                   _elemtype(et)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn = E*h_/((1-nu)*(1+nu));
  }
  virtual ~IsotropicElasticForceBulkTermC0Plate(){}
  virtual void get(MElement *ele,int npts,IntPt *GP,fullVector<double> &v);
}; // IsotropicElasticStiffBulkTermC0Plate


class IsotropicElasticStiffBulkTermC0Plate : public DgC0BilinearTerm<SVector3,SVector3>
{
 protected :
  double E,nu,h,Cm,Cn;
  bool sym,fullDg;
  displacementField *ufield;
  IPField *ipf;
  DgC0LinearTerm<SVector3> *lterm;

 public :
  IsotropicElasticStiffBulkTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,DgC0FunctionSpace<SVector3>& space2_,materialLaw* mlaw,
                                  double h_,bool FullDG, displacementField *uf, IPField *ip, SolElementType::eltype elt) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space2_),h(h_),
                                                            fullDg(FullDG), ufield(uf), ipf(ip)
  {
    sym=(&space1_==&space2_);
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn = E*h_/((1-nu)*(1+nu));
    lterm = new IsotropicElasticForceBulkTermC0Plate(space1,mlaw,h,fullDg,ufield,ipf,elt);
  }
  IsotropicElasticStiffBulkTermC0Plate(DgC0FunctionSpace<SVector3>& space1_, materialLaw* mlaw,
                                  double h_,bool FullDG, displacementField *uf,IPField *ip, SolElementType::eltype elt) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space1_),
                                                            h(h_),fullDg(FullDG), ufield(uf), ipf(ip)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    sym=true;
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn = E*h_/((1-nu)*(1+nu));
    lterm = new IsotropicElasticForceBulkTermC0Plate(space1,mlaw,h,fullDg,ufield,ipf,elt);
  }
  virtual ~IsotropicElasticStiffBulkTermC0Plate(){}
  virtual void get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m);
  DgC0LinearTerm<SVector3>* getLinearTerm() const{return lterm;}
}; // IsotropicElasticStiffBulkTermC0Plate

class IsotropicElasticForceInterfaceTermC0Plate : public DgC0LinearTerm<SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool  fullDg;
  displacementField *ufield;
  IPField *ipf;
  SolElementType::eltype _elemtype; // TODO Two SolElementType One for minus and one for + element
  bool MatrixByPerturbation, minus; // Use to compute matrix by numerical perturbation
  Dof pertDof;
  double pert;
 public :
  IsotropicElasticForceInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_, materialLaw *mlaw, double beta1_,
                                       double beta2_, double beta3_, double h_, bool FullDg_, displacementField *uf,
                                       IPField *ip,
                                       SolElementType::eltype et,bool mbp=false) : DgC0LinearTerm<SVector3>(space1_),
                                                                            beta1(beta1_), beta2(beta2_),beta3(beta3_),h(h_),
                                                                            fullDg(FullDg_), ufield(uf), ipf(ip),
                                                                            _elemtype(et), MatrixByPerturbation(mbp), minus(true),
                                                                            pert(0.), pertDof(0,0)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn =  E*h_/((1-nu)*(1+nu));
    Cs = ( E * h_ ) / (2 * ( 1 + nu ) );
  }
  virtual ~IsotropicElasticForceInterfaceTermC0Plate() {}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullVector<double> &v);
  void setMinus(const bool e){minus = e;}
  virtual void setPert(const double eps){pert = eps;}
  virtual void setPertDof(const Dof &D){pertDof = Dof(D.getEntity(),D.getType());}
  virtual void invSign(){inverseSign ? inverseSign = false : inverseSign=true;}
}; // class IsotropicElasticForceInterfaceTermC0Plate

class IsotropicElasticStiffInterfaceTermC0Plate : public DgC0BilinearTerm<SVector3,SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool sym;
  bool  fullDg;
  SolElementType::eltype _elemtype; // TODO 2 SolElementType (one for minus and one for plus element)
  displacementField *ufield;
  IPField *ipf;
  double perturbation; // value of perturbation for compution by numerical perturbation.
  DgC0LinearTerm<SVector3> *lterm;

 public :
  IsotropicElasticStiffInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,DgC0FunctionSpace<SVector3>& space2_,materialLaw* mlaw,
                                       double beta1_,double beta2_, double beta3_, double h_, displacementField *uf,
                                       IPField *ip, SolElementType::eltype elt_,
                                       bool FullDg_, double eps=1.e-8) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space2_),beta1(beta1_),
                                                        beta2(beta2_),beta3(beta3_),h(h_), fullDg(FullDg_), _elemtype(elt_), ufield(uf),
                                                        ipf(ip), perturbation(eps)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn =  E*h_/((1-nu)*(1+nu));
    Cs = ( E * h_ ) / (2 * ( 1 + nu ) );
    sym=(&space1_==&space2_);
    lterm = new IsotropicElasticForceInterfaceTermC0Plate (space1,mlaw,beta1,beta2,beta3,h,fullDg,ufield,ipf,_elemtype,true);
  }
  IsotropicElasticStiffInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,materialLaw *mlaw,
                                       double beta1_,double beta2_, double beta3_, double h_, displacementField *uf,
                                       IPField *ip, SolElementType::eltype elt_,
                                       bool FullDg_, double eps=1.e-8) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space1_),
                                                     beta1(beta1_),beta2(beta2_),beta3(beta3_),h(h_),fullDg(FullDg_), _elemtype(elt_),
                                                     ufield(uf), ipf(ip), perturbation(eps)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn =  E*h_/((1-nu)*(1+nu));
    Cs = ( E * h_ ) / (2 * ( 1 + nu ) );
    sym=true;
    lterm = new IsotropicElasticForceInterfaceTermC0Plate(space1,mlaw,beta1,beta2,beta3,h,fullDg,ufield,ipf,_elemtype,true);
  }
  virtual ~IsotropicElasticStiffInterfaceTermC0Plate() {delete lterm;}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullMatrix<double> &m);
  IsotropicElasticForceInterfaceTermC0Plate* getLinearTerm() {return dynamic_cast<IsotropicElasticForceInterfaceTermC0Plate*>(lterm);}
}; // class IsotropicElasticStiffInterfaceTermC0Plate


class IsotropicElasticForceVirtualInterfaceTermC0Plate : public DgC0LinearTerm<SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool  fullDg;
  displacementField *ufield;
  IPField *ipf;
  SolElementType::eltype _elemtype;

 public :
  IsotropicElasticForceVirtualInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,materialLaw *mlaw, double beta1_,
                                       double beta2_, double beta3_, double h_, bool FullDg_, displacementField *uf,
                                       IPField  *ip,
                                       SolElementType::eltype et) : DgC0LinearTerm<SVector3>(space1_),
                                                                                            beta1(beta1_),
                                                                                            beta2(beta2_),beta3(beta3_),h(h_),
                                                                                            fullDg(FullDg_),
                                                                                            ufield(uf),ipf(ip), _elemtype(et)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn =  E*h_/((1-nu)*(1+nu));
    Cs = ( E * h_ ) / (2 * ( 1 + nu ) );
  }
  virtual ~IsotropicElasticForceVirtualInterfaceTermC0Plate() {}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullVector<double> &v);
  virtual void invSign(){inverseSign ? inverseSign = false : inverseSign=true;}
}; // class IsotropicElasticForceVirtualInterfaceTermC0Plate


class IsotropicElasticStiffVirtualInterfaceTermC0Plate : public DgC0BilinearTerm<SVector3,SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool sym;
  bool  fullDg;
  displacementField *ufield;
  IPField *ipf;
  DgC0LinearTerm<SVector3> *lterm;
  SolElementType::eltype et;

 public :
  IsotropicElasticStiffVirtualInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,DgC0FunctionSpace<SVector3>& space2_,materialLaw *mlaw,
                                       double beta1_,double beta2_, double beta3_, double h_, displacementField *uf, IPField *ipf_, SolElementType::eltype et_,
                                       bool FullDg_) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space2_),beta1(beta1_),
                                                        beta2(beta2_),beta3(beta3_),h(h_), fullDg(FullDg_) ,ufield(uf),ipf(ipf_), et(et_)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn =  E*h_/((1-nu)*(1+nu));
    Cs = ( E * h_ ) / (2 * ( 1 + nu ) );
    sym=(&space1_==&space2_);
    lterm = new IsotropicElasticForceVirtualInterfaceTermC0Plate(space1,mlaw,beta1,beta2,beta3,h,fullDg,ufield,ipf,et);
  }
  IsotropicElasticStiffVirtualInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,materialLaw* mlaw,
                                       double beta1_,double beta2_, double beta3_,
                                       double h_, displacementField *uf, IPField *ipf_,SolElementType::eltype et_,
                                       bool FullDg_) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space1_),
                                                     beta1(beta1_),beta2(beta2_),beta3(beta3_),h(h_),fullDg(FullDg_) ,ufield(uf), ipf(ipf_), et(et_)
  {
    linearElasticLawPlaneStress *mlt = dynamic_cast<linearElasticLawPlaneStress*>(mlaw);
    E = mlt->getYoung();
    nu = mlt->getPoisson();
    Cm = ( E * h_ * h_ * h_ ) / ( 12 * (1 - nu) * (1 + nu) );
    Cn =  E*h_/((1-nu)*(1+nu));
    Cs = ( E * h_ ) / (2 * ( 1 + nu ) );
    sym=true;
    lterm = new IsotropicElasticForceVirtualInterfaceTermC0Plate(space1,mlaw,beta1,beta2,beta3,h,fullDg,ufield,ipf,et);
  }
  virtual ~IsotropicElasticStiffVirtualInterfaceTermC0Plate() {delete lterm;}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullMatrix<double> &m);
  DgC0LinearTerm<SVector3>* getLinearTerm() const{return lterm;}
}; // class IsotropicElasticStiffVirtualInterfaceTermC0Plate

#endif // DGC0PLATETERMS_H
