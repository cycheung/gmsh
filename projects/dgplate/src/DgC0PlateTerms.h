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

class  DgC0BilinearTermBase
{
 public :
  virtual ~DgC0BilinearTermBase() {}
  virtual void get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m) = 0 ;
  //virtual void get(MElement *ele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m) = 0 ;
  virtual void getKeys(MElement*, std::vector<Dof>&)=0; // To have the Dof with different formulation
/*  virtual void getInter(MInterfaceElement *iele,int npts,IntPt *GP,bool FullDg, fullMatrix<double> &m){} ;
  virtual void getInterForce(MInterfaceElement *iele,int npts,IntPt *GP,const fullMatrix<double> &disp,fullMatrix<double> &m){} ;
  virtual void getForce(MElement *ele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m){}
  virtual void getInterOnBoundary(MInterfaceElement *iele,int npts,IntPt *GP,fullMatrix<double> &m){}
  virtual void getInterForceOnBoundary(MInterfaceElement *iele,int npts,IntPt *GP,const fullMatrix<double> &disp,fullMatrix<double> &m){} ;*/
};

template<class T1,class T2> class DgC0BilinearTerm : public DgC0BilinearTermBase
{
 protected :
  DgC0FunctionSpace<T1>& space1;
  DgC0FunctionSpace<T2>& space2;
 public :
  DgC0BilinearTerm(DgC0FunctionSpace<T1>& space1_,DgC0FunctionSpace<T1>& space2_) : space1(space1_),space2(space2_) {}
  virtual ~DgC0BilinearTerm() {}
};

class  DgC0LinearTermBase
{
  public:
  virtual ~DgC0LinearTermBase() {}
  virtual void get(MElement *ele,int npts,IntPt *GP,fullVector<double> &v) =0;
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys)=0;
};

template<class T1> class DgC0LinearTerm : public DgC0LinearTermBase
{
 protected :
  DgC0FunctionSpace<T1>& space1;
 public :
  DgC0LinearTerm(DgC0FunctionSpace<T1>& space1_) : space1(space1_) {}
  virtual ~DgC0LinearTerm() {}
};

class  DgC0ScalarTermBase
{
 public :
  virtual ~DgC0ScalarTermBase() {}
  virtual void get(MElement *ele,int npts,IntPt *GP,double &val) =0;
};

class DgC0ScalarTerm : public DgC0ScalarTermBase
{
 public :
  virtual ~DgC0ScalarTerm() {}
};


template<class T1,class T2> class DgC0BilinearTermToScalarTerm : public DgC0ScalarTerm
{
  DgC0BilinearTerm<T1,T2> &bilterm;
  public :
  DgC0BilinearTermToScalarTerm(DgC0BilinearTerm<T1,T2> &bilterm_): bilterm(bilterm_){}
  virtual ~DgC0BilinearTermToScalarTerm() {}
  virtual void get(MElement *ele,int npts,IntPt *GP,double &val)
  {
    fullMatrix<double> localMatrix;
    bilterm.get(ele,npts,GP,localMatrix);
    val=localMatrix(0,0);
  }
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){bilterm.getKeys(ele,keys);}
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
  virtual void getKeys(MElement *ele,std::vector<Dof> &keys){DgC0LinearTerm<T1>::space1.getKeys(ele,keys);}
};

class IsotropicElasticForceBulkTermC0Plate : public DgC0LinearTerm<SVector3>
{
 protected :
  double E,nu,h,Cm,Cn;
  bool fullDg,inverseSign;
  displacementField *ufield;
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf;
  SolElementType::eltype _elemtype;

 public :
  IsotropicElasticForceBulkTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,double E_,double nu_,double h_, bool FullDG,
                                       displacementField *uf,IPField<DGelasticField, DgC0FunctionSpace<SVector3> >*ip,
                                       SolElementType::eltype et,bool no=false) : DgC0LinearTerm<SVector3>(space1_),E(E_),nu(nu_),h(h_),
                                                                                   fullDg(FullDG), inverseSign(no), ufield(uf), ipf(ip),
                                                                                   _elemtype(et)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn = E_*h_/((1-nu_)*(1+nu_));
  }
  virtual ~IsotropicElasticForceBulkTermC0Plate(){}
  virtual void get(MElement *ele,int npts,IntPt *GP,fullVector<double> &v);
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){space1.getKeys(ele,keys,fullDg);}
}; // IsotropicElasticStiffBulkTermC0Plate

class IsotropicElasticStiffBulkTermC0Plate : public DgC0BilinearTerm<SVector3,SVector3>
{
 protected :
  double E,nu,h,Cm,Cn;
  bool sym,fullDg;
  //displacementField *ufield;
  //DgC0LinearTerm<SVector3> *lterm;

 public :
  IsotropicElasticStiffBulkTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,DgC0FunctionSpace<SVector3>& space2_,double E_,double nu_,
                                  double h_,bool FullDG) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space2_),E(E_),nu(nu_),h(h_),
                                                            fullDg(FullDG)
  {
    sym=(&space1_==&space2_);
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn = E_*h_/((1-nu_)*(1+nu_));
    //lterm = new IsotropicElasticForceBulkTermC0Plate(space1,E,nu,h,fullDg,ufield);
  }
  IsotropicElasticStiffBulkTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,double E_,double nu_,
                                  double h_,bool FullDG) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space1_),
                                                            E(E_),nu(nu_),h(h_),fullDg(FullDG)
  {
    sym=true;
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn = E_*h_/((1-nu_)*(1+nu_));
    //lterm = new IsotropicElasticForceBulkTermC0Plate(space1,E,nu,h,fullDg,ufield);
  }
  virtual ~IsotropicElasticStiffBulkTermC0Plate(){}
  virtual void get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m);
  //virtual void get(MElement *ele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m);
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){space1.getKeys(ele,keys,fullDg);}
  //DgC0LinearTerm<SVector3>* getLinearTerm() const{return lterm;}
}; // IsotropicElasticStiffBulkTermC0Plate

class IsotropicElasticForceInterfaceTermC0Plate : public DgC0LinearTerm<SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool  fullDg, inverseSign;
  displacementField *ufield;
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf;
  SolElementType::eltype _elemtype; // TODO Two SolElementType One for minus and one for + element

 public :
  IsotropicElasticForceInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_, double E_,double nu_, double beta1_,
                                       double beta2_, double beta3_, double h_, bool FullDg_, displacementField *uf,
                                       IPField<DGelasticField, DgC0FunctionSpace<SVector3> >  *ip,
                                       SolElementType::eltype et, bool no=false) : DgC0LinearTerm<SVector3>(space1_),E(E_),nu(nu_),
                                                                            beta1(beta1_), beta2(beta2_),beta3(beta3_),h(h_),
                                                                            fullDg(FullDg_), inverseSign(no), ufield(uf), ipf(ip),
                                                                            _elemtype(et)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn =  E_*h_/((1-nu_)*(1+nu_));
    Cs = ( E_ * h_ ) / (2 * ( 1 + nu_ ) );
  }
  virtual ~IsotropicElasticForceInterfaceTermC0Plate() {}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullVector<double> &v);
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){
    MInterfaceElement *ie = dynamic_cast<MInterfaceElement*>(ele);
    space1.getKeys(ie->getElem(0),keys,fullDg);
    space1.getKeys(ie->getElem(1),keys,fullDg);
  }

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
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf;
  double perturbation; // value of perturbation for compution by numerical perturbation.
  DgC0LinearTerm<SVector3> *lterm;

 public :
  IsotropicElasticStiffInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,DgC0FunctionSpace<SVector3>& space2_,double E_,double nu_,
                                       double beta1_,double beta2_, double beta3_, double h_, displacementField *uf,
                                       IPField<DGelasticField, DgC0FunctionSpace<SVector3> >  *ip, SolElementType::eltype elt_,
                                       bool FullDg_, double eps=1.e-8) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space2_),E(E_),nu(nu_),beta1(beta1_),
                                                        beta2(beta2_),beta3(beta3_),h(h_), fullDg(FullDg_), _elemtype(elt_), ufield(uf),
                                                        ipf(ip), perturbation(eps)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn =  E_*h_/((1-nu_)*(1+nu_));
    Cs = ( E_ * h_ ) / (2 * ( 1 + nu_ ) );
    sym=(&space1_==&space2_);
    lterm = new IsotropicElasticForceInterfaceTermC0Plate (space1,E,nu,beta1,beta2,beta3,h,fullDg,ufield,ipf,_elemtype,false);
  }
  IsotropicElasticStiffInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,double E_,double nu_,
                                       double beta1_,double beta2_, double beta3_, double h_, displacementField *uf,
                                       IPField<DGelasticField, DgC0FunctionSpace<SVector3> >  *ip, SolElementType::eltype elt_,
                                       bool FullDg_, double eps=1.e-8) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space1_),E(E_),nu(nu_),
                                                     beta1(beta1_),beta2(beta2_),beta3(beta3_),h(h_),fullDg(FullDg_), _elemtype(elt_),
                                                     ufield(uf), ipf(ip), perturbation(eps)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn =  E_*h_/((1-nu_)*(1+nu_));
    Cs = ( E_ * h_ ) / (2 * ( 1 + nu_ ) );
    sym=true;
    lterm = new IsotropicElasticForceInterfaceTermC0Plate(space1,E,nu,beta1,beta2,beta3,h,fullDg,ufield,ipf,_elemtype,false);
  }
  virtual ~IsotropicElasticStiffInterfaceTermC0Plate() {}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullMatrix<double> &m);
  //virtual void get(MElement *iele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m);
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){
    MInterfaceElement *ie = dynamic_cast<MInterfaceElement*>(ele);
    space1.getKeys(ie->getElem(0),keys,fullDg);
    space1.getKeys(ie->getElem(1),keys,fullDg);
  }
  DgC0LinearTerm<SVector3>* getLinearTerm() const{return lterm;}
}; // class IsotropicElasticStiffInterfaceTermC0Plate


class IsotropicElasticForceVirtualInterfaceTermC0Plate : public DgC0LinearTerm<SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool  fullDg, inverseSign;
  displacementField *ufield;
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> >  *ipf;
  SolElementType::eltype _elemtype;

 public :
  IsotropicElasticForceVirtualInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,double E_,double nu_, double beta1_,
                                       double beta2_, double beta3_, double h_, bool FullDg_, displacementField *uf,
                                       IPField<DGelasticField,DgC0FunctionSpace<SVector3> >  *ip,
                                       SolElementType::eltype et, bool no=false) : DgC0LinearTerm<SVector3>(space1_),
                                                                                            E(E_),nu(nu_),beta1(beta1_),
                                                                                            beta2(beta2_),beta3(beta3_),h(h_),
                                                                                            fullDg(FullDg_), inverseSign(no),
                                                                                            ufield(uf),ipf(ip), _elemtype(et)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn =  E_*h_/((1-nu_)*(1+nu_));
    Cs = ( E_ * h_ ) / (2 * ( 1 + nu_ ) );
  }
  virtual ~IsotropicElasticForceVirtualInterfaceTermC0Plate() {}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullVector<double> &v);
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){
    MInterfaceElement *ie = dynamic_cast<MInterfaceElement*>(ele);
    space1.getKeys(ie->getElem(0),keys,fullDg);
  }
}; // class IsotropicElasticForceVirtualInterfaceTermC0Plate


class IsotropicElasticStiffVirtualInterfaceTermC0Plate : public DgC0BilinearTerm<SVector3,SVector3>
{
 protected :
  double E,nu,beta1,beta2,beta3,h;
  double Cm,Cn,Cs;
  bool sym;
  bool  fullDg;
  displacementField *ufield;
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf;
  DgC0LinearTerm<SVector3> *lterm;
  bool no;
  SolElementType::eltype et;

 public :
  IsotropicElasticStiffVirtualInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,DgC0FunctionSpace<SVector3>& space2_,double E_,double nu_,
                                       double beta1_,double beta2_, double beta3_, double h_, displacementField *uf, IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf_, SolElementType::eltype et_, bool no_,
                                       bool FullDg_) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space2_),E(E_),nu(nu_),beta1(beta1_),
                                                        beta2(beta2_),beta3(beta3_),h(h_), fullDg(FullDg_) ,ufield(uf),ipf(ipf_), et(et_), no(no_)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn =  E_*h_/((1-nu_)*(1+nu_));
    Cs = ( E_ * h_ ) / (2 * ( 1 + nu_ ) );
    sym=(&space1_==&space2_);
    lterm = new IsotropicElasticForceVirtualInterfaceTermC0Plate(space1,E,nu,beta1,beta2,beta3,h,fullDg,ufield,ipf,et,false);
  }
  IsotropicElasticStiffVirtualInterfaceTermC0Plate(DgC0FunctionSpace<SVector3>& space1_,double E_,double nu_,
                                       double beta1_,double beta2_, double beta3_,
                                       double h_, displacementField *uf, IPField<DGelasticField,DgC0FunctionSpace<SVector3> > *ipf_,SolElementType::eltype et_, bool no_, bool FullDg_) : DgC0BilinearTerm<SVector3,SVector3>(space1_,space1_),E(E_),nu(nu_),
                                                     beta1(beta1_),beta2(beta2_),beta3(beta3_),h(h_),fullDg(FullDg_) ,ufield(uf), ipf(ipf_), et(et_)
  {
    Cm = ( E_ * h_ * h_ * h_ ) / ( 12 * (1 - nu_) * (1 + nu_) );
    Cn =  E_*h_/((1-nu_)*(1+nu_));
    Cs = ( E_ * h_ ) / (2 * ( 1 + nu_ ) );
    sym=true;
    lterm = new IsotropicElasticForceVirtualInterfaceTermC0Plate(space1,E,nu,beta1,beta2,beta3,h,fullDg,ufield,ipf,et,false);
  }
  virtual ~IsotropicElasticStiffVirtualInterfaceTermC0Plate() {}
  virtual void get(MElement *iele,int npts,IntPt *GP, fullMatrix<double> &m);
  //virtual void get(MElement *iele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m);
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){
    MInterfaceElement *ie = dynamic_cast<MInterfaceElement*>(ele);
    space1.getKeys(ie->getElem(0),keys,fullDg);
  }
  DgC0LinearTerm<SVector3>* getLinearTerm() const{return lterm;}
}; // class IsotropicElasticStiffVirtualInterfaceTermC0Plate

#endif // DGC0PLATETERMS_H
