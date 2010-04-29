//
// C++ Interface: terms
//
// Description: FunctionSpace used (scalar lagragian function space in 3D)
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef DGC0PLATEFUNCTIONSPACE_H_
#define DGC0PLATEFUNCTIONSPACE_H_
#include "MElement.h"
#include "functionSpace.h"
#include "MInterfaceElement.h"
#include "DgC0PlateDof.h"

class DgC0FunctionSpaceBase
{
 public:
  virtual int getNumKeys(MElement *ele)=0; // if one needs the number of dofs
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){printf("Warning not numbered dof");} //not pure virtual because this one is not definite
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys, bool FullDg, std::vector<int> *numvertex=NULL) {printf("Warning not numbered dof");}  // for the space where this one is used
};

template<class T>
class DgC0FunctionSpace : public DgC0FunctionSpaceBase
{
 public:
  typedef typename TensorialTraits<double>::ValType ValType;
  typedef typename TensorialTraits<double>::GradType GradType;
  typedef typename TensorialTraits<double>::HessType HessType;
  virtual void f(MElement *ele, double u, double v, double w, std::vector<ValType> &vals){};
  virtual void gradf(MElement *ele, double u, double v, double w,std::vector<GradType> &grads){};
  virtual void fuvw(MElement *ele, double u, double v, double w, std::vector<ValType> &vals){};
  virtual void gradfuvw(MElement *ele, double u, double v, double w,std::vector<GradType> &grads) {}; // should return to pure virtual once all is done.
  virtual void hessfuvw(MElement *ele, double u, double v, double w,std::vector<HessType> &hess)=0;
  virtual int getNumKeys(MElement *ele)=0; // if one needs the number of dofs
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){printf("Warning not implemented function\n");};
  virtual int getId(void) const{printf("Warning not implemented function\n");return 0;}
};

class DgC0ScalarLagrangeFunctionSpace : public DgC0FunctionSpace<double>
{
 public:
  typedef TensorialTraits<double>::ValType ValType;
  typedef TensorialTraits<double>::GradType GradType;
  typedef TensorialTraits<double>::HessType HessType;

 protected:
  int _iField; // field number (used to build dof keys)

 private:
  virtual void getKeys(MVertex *ver, std::vector<Dof> &keys)
  {
    keys.push_back(Dof(ver->getNum(), _iField));
  }
  virtual void getKeys(const int numelem,const int vernum, std::vector<Dof> &keys)
  {
    keys.push_back(Dof(numelem, DgC0PlateDof::createTypeWithThreeInts(0,_iField,vernum)));
  }

 public:
  DgC0ScalarLagrangeFunctionSpace(int i=0):_iField(i) {}
  int getId(void) const {return _iField;};
  virtual void f(MElement *ele, double u, double v, double w, std::vector<ValType> &vals)
  {
    if (ele->getParent()) ele = ele->getParent();
    int ndofs= ele->getNumVertices();
    int curpos=vals.size();
    vals.resize(curpos+ndofs);
    ele->getShapeFunctions(u, v, w, &(vals[curpos]));
  }


  // Fonction renvoyant un vecteur contenant le grandient de chaque FF
  virtual void gradf(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)
  {
    if (ele->getParent()) ele = ele->getParent();
    int ndofs= ele->getNumVertices();
    grads.reserve(grads.size()+ndofs);
    double gradsuvw[256][3];
    ele->getGradShapeFunctions(u, v, w, gradsuvw);
    double jac[3][3];
    double invjac[3][3];
    const double detJ = ele->getJacobian(u, v, w, jac); // redondant : on fait cet appel a l'exterieur
    inv3x3(jac, invjac);
    for (int i=0;i<ndofs;++i)
      grads.push_back(GradType(
      invjac[0][0] * gradsuvw[i][0] + invjac[0][1] * gradsuvw[i][1] + invjac[0][2] * gradsuvw[i][2],
      invjac[1][0] * gradsuvw[i][0] + invjac[1][1] * gradsuvw[i][1] + invjac[1][2] * gradsuvw[i][2],
      invjac[2][0] * gradsuvw[i][0] + invjac[2][1] * gradsuvw[i][1] + invjac[2][2] * gradsuvw[i][2]
                            ));
  };
    // Fonction renvoyant un vecteur contenant le hessien [][] de chaque FF dans l'espace ISOPARAMETRIQUE
  virtual void hessfuvw(MElement *ele, double u, double v, double w,std::vector<HessType> &hess)
  {
    int ndofs= ele->getNumVertices();  // ATTENTION RETOURNE LE NBBRE DE NOEUDS ET PAS LE NBRE DE DDL
    hess.reserve(hess.size()+ndofs);   // permet de mettre les composantes suivantes à la suite du vecteur
    double hessuvw[256][3][3];
    ele->getHessShapeFunctions(u, v, w, hessuvw);
    HessType hesst;
    for (int i=0;i<ndofs;++i){
      hesst(0,0) = hessuvw[i][0][0] ; hesst(0,1) = hessuvw[i][0][1] ; hesst(0,2) = hessuvw[i][0][2];
      hesst(1,0) = hessuvw[i][1][0] ; hesst(1,1) = hessuvw[i][1][1] ; hesst(1,2) = hessuvw[i][1][2];
      hesst(2,0) = hessuvw[i][2][0] ; hesst(2,1) = hessuvw[i][2][1] ; hesst(2,2) = hessuvw[i][2][2];
      hess.push_back(hesst);
    }
  };

  virtual void gradfuvw(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)
  {
    if (ele->getParent()) ele = ele->getParent();
    int ndofs= ele->getNumVertices();
    grads.reserve(grads.size()+ndofs);
    double gradsuvw[256][3];
    ele->getGradShapeFunctions(u, v, w, gradsuvw);
    for (int i=0;i<ndofs;++i)
      grads.push_back(GradType(
      gradsuvw[i][0] ,
      gradsuvw[i][1] ,
      gradsuvw[i][2]
                      ));
  }

  virtual void fuvw(MElement *ele, double u, double v, double w,std::vector<ValType> &vals)
  {
    if (ele->getParent()) ele = ele->getParent();
    int ndofs= ele->getNumVertices();
    vals.reserve(vals.size()+ndofs);
    double valsuvw[256];
    ele->getShapeFunctions(u, v, w, valsuvw);
    for (int i=0;i<ndofs;++i) {vals.push_back(valsuvw[i]);}
  }

  virtual int getNumKeys(MElement *ele) {if (ele->getParent()) ele = ele->getParent();return ele->getNumVertices();}

  virtual void getKeys(MElement *ele, std::vector<Dof> &keys) // appends ...
  {
      if (ele->getParent()) ele = ele->getParent();
      int ndofs= ele->getNumVertices();
      keys.reserve(keys.size()+ndofs);
      for (int i=0;i<ndofs;++i)
        getKeys(ele->getVertex(i),keys);
  }
    void getKeys(MElement *ele, std::vector<Dof> &keys,bool FullDg,int vernum=0) // appends ... vernum used for FullDg only
  {
      if (ele->getParent()) ele = ele->getParent();
      int ndofs= ele->getNumVertices();
      keys.reserve(keys.size()+ndofs);
      if(!FullDg){
        for (int i=0;i<ndofs;++i)
          getKeys(ele->getVertex(i),keys);
      }
      else{
        for (int i=0;i<ndofs;++i)
          getKeys(ele->getNum(),vernum,keys);
      }
  }
};

template <class T> class DgC0ScalarToAnyFunctionSpace : public DgC0FunctionSpace<T>
{
protected :
  std::vector<int> comp;
public :
  typedef typename TensorialTraits<double>::ValType ValType;
  typedef typename TensorialTraits<double>::GradType GradType;
  typedef typename TensorialTraits<double>::HessType HessType;
protected :
  DgC0FunctionSpace<double> *ScalarFS;
public :
  template <class T2> DgC0ScalarToAnyFunctionSpace(const T2 &SFS, int comp1_): ScalarFS(new T2(SFS))
  {comp.push_back(comp1_);}
  template <class T2> DgC0ScalarToAnyFunctionSpace(const T2 &SFS, int comp1_, int comp2_): ScalarFS(new T2(SFS))
  {comp.push_back(comp1_),comp.push_back(comp2_);}
  template <class T2> DgC0ScalarToAnyFunctionSpace(const T2 &SFS, int comp1_, int comp2_, int comp3_): ScalarFS(new T2(SFS))
  {comp.push_back(comp1_),comp.push_back(comp2_),comp.push_back(comp3_);}

  virtual ~DgC0ScalarToAnyFunctionSpace() {delete ScalarFS;}

  virtual void f(MElement *ele, double u, double v, double w, std::vector<ValType> &vals)
  {
    ScalarFS->f(ele,u,v,w,vals);
  }

  virtual void gradf(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)
  {
    ScalarFS->gradf(ele,u,v,w,grads);
  }
  virtual void hessfuvw(MElement *ele, double u, double v, double w,std::vector<HessType> &hess)
  {
    ScalarFS->hessfuvw(ele,u,v,w,hess);
  }
  virtual void gradfuvw(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)
  {
    ScalarFS->gradfuvw(ele,u,v,w,grads);
  }

    virtual void fuvw(MElement *ele, double u, double v, double w,std::vector<ValType> &vals)
  {
    ScalarFS->fuvw(ele,u,v,w,vals);
  }

  virtual int getNumKeys(MElement *ele) {return ScalarFS->getNumKeys(ele)*comp.size();}

  virtual void getKeys(MElement *ele, std::vector<Dof> &keys, bool FullDg,std::vector<int> *numvertex=NULL)
  {
    int nk=ScalarFS->getNumKeys(ele); // return the number of vertices
    std::vector<Dof> bufk;
    bufk.reserve(nk);
    ScalarFS->getKeys(ele,bufk);
    int nbdofs=bufk.size();
    int nbcomp=comp.size();
    int curpos=keys.size();
    keys.reserve(curpos+nbcomp*nbdofs);
    if(!FullDg){
      for (int j=0;j<nbcomp;++j)
      {
        for (int i=0;i<nk;++i)
        {
          int i1,i2;
          DgC0PlateDof::getThreeIntsFromType(bufk[i].getType(), i1,i2);
          keys.push_back(Dof(bufk[i].getEntity(),DgC0PlateDof::createTypeWithThreeInts(comp[j],i1)));
        }
      }
    }
    else{
      if(numvertex==NULL){ // All vertex are generated
        for (int j=0;j<nbcomp;++j)
        {
          for (int i=0;i<nk;i++)
          {
            int i1,i2,i3;
            DgC0PlateDof::getThreeIntsFromType(bufk[i].getType(), i1,i2);
            keys.push_back(Dof(ele->getNum(),DgC0PlateDof::createTypeWithThreeInts(comp[j],i1,i)));
          }
        }
      }
      else{ // Only vertex given in numvertex are generated
        for (int j=0;j<nbcomp;++j)
        {
          for (int i=0;i<(*numvertex).size();i++)
          {
            int i1,i2,i3;
            DgC0PlateDof::getThreeIntsFromType(bufk[i].getType(), i1,i2);
            keys.push_back(Dof(ele->getNum(),DgC0PlateDof::createTypeWithThreeInts(comp[j],i1,(*numvertex)[i])));
          }
        }
      }
    }
  }
};

class DgC0LagrangeFunctionSpace : public DgC0ScalarToAnyFunctionSpace<SVector3>
{
 public:
  enum Along { VECTOR_X=0, VECTOR_Y=1, VECTOR_Z=2 };
  typedef TensorialTraits<double>::ValType ValType;
  typedef TensorialTraits<double>::GradType GradType;
  typedef TensorialTraits<double>::HessType HessType;
  DgC0LagrangeFunctionSpace(int id) :
          DgC0ScalarToAnyFunctionSpace<SVector3>::DgC0ScalarToAnyFunctionSpace(DgC0ScalarLagrangeFunctionSpace(id),0,1,2)
  {}
  DgC0LagrangeFunctionSpace(int id,Along comp1) :
          DgC0ScalarToAnyFunctionSpace<SVector3>::DgC0ScalarToAnyFunctionSpace(DgC0ScalarLagrangeFunctionSpace(id),comp1)
  {}
  DgC0LagrangeFunctionSpace(int id,Along comp1,Along comp2) :
          DgC0ScalarToAnyFunctionSpace<SVector3>::DgC0ScalarToAnyFunctionSpace(DgC0ScalarLagrangeFunctionSpace(id),comp1,comp2)
  {}
  DgC0LagrangeFunctionSpace(int id,Along comp1,Along comp2, Along comp3) :
          DgC0ScalarToAnyFunctionSpace<SVector3>::DgC0ScalarToAnyFunctionSpace(DgC0ScalarLagrangeFunctionSpace(id),comp1,comp2,comp3)
  {}
};

#endif // DGC0PLATEFUNCTIONSPACE_H_
