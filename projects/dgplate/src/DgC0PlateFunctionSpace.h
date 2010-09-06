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

class DgC0FunctionSpaceBase : public FunctionSpaceBase{
 public :
  virtual void getKeys(MInterfaceElement *iele, std::vector<Dof> &keys)=0;
};

class DgC0FullDgFunctionSpaceBase
{
 public:
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys, std::vector<int> *numvertex)=0;
};

template<class T>
class DgC0FunctionSpace : public DgC0FunctionSpaceBase,public FunctionSpace<T>
{
 public:
  typedef TensorialTraits<double>::ValType ValType;
  typedef TensorialTraits<double>::GradType GradType;
  typedef TensorialTraits<double>::HessType HessType;
  virtual void f(MElement *ele, double u, double v, double w, std::vector<ValType> &vals)=0;
  virtual void gradf(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)=0;
  virtual void fuvw(MElement *ele, double u, double v, double w, std::vector<ValType> &vals)=0;
  virtual void gradfuvw(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)=0;
  virtual void hessfuvw(MElement *ele, double u, double v, double w,std::vector<HessType> &hess)=0;
  virtual int getNumKeys(MElement *ele)=0; // if one needs the number of dofs
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys)=0;
  virtual int getId(void) const=0;//{printf("Warning not implemented function\n");return 0;}
  // Same definition for all Function Space (search the keys of two elements)
  virtual void getKeys(MInterfaceElement *iele, std::vector<Dof> &keys){
    this->getKeys(iele->getElem(0),keys);
    if(!(iele->getElem(0) == iele->getElem(1)))
      this->getKeys(iele->getElem(1),keys);
  }
};

class DgC0LagrangeFunctionSpace : public DgC0FunctionSpace<SVector3>
{
 protected :
//  int iField; // otherwise getId in scalarLagrangeFunctionSpace doesn't work ?? //TODO How to access
  std::vector<int> comp;
  ScalarLagrangeFunctionSpace *ScalarFS;
 public:
  typedef TensorialTraits<double>::ValType ValType;
  typedef TensorialTraits<double>::GradType GradType;
  typedef TensorialTraits<double>::HessType HessType;

  DgC0LagrangeFunctionSpace(int id) : ScalarFS(new ScalarLagrangeFunctionSpace(id))
  {comp.push_back(0),comp.push_back(1),comp.push_back(2);}
  ~DgC0LagrangeFunctionSpace(){delete ScalarFS;}
  virtual void f(MElement *ele, double u, double v, double w, std::vector<SVector3> &vals){ };
  virtual void gradf(MElement *ele, double u, double v, double w,std::vector<STensor3> &grads){ };

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

  virtual int getNumKeys(MElement *ele) {if (ele->getParent()) ele = ele->getParent();return comp.size()*ele->getNumVertices();}
  virtual int getId(void) const {return ScalarFS->getId();}

  virtual void getKeys(MElement *ele, std::vector<Dof> &keys)=0;
};

class DgC0CgDgLagrangeFunctionSpace : public DgC0LagrangeFunctionSpace{
 public :
  DgC0CgDgLagrangeFunctionSpace(int id) : DgC0LagrangeFunctionSpace(id){};
  ~DgC0CgDgLagrangeFunctionSpace(){};

  virtual void getKeys(MElement *ele, std::vector<Dof> &keys)
  {
    int nk=ScalarFS->getNumKeys(ele); // return the number of vertices
    std::vector<Dof> bufk;
    bufk.reserve(nk);
    ScalarFS->getKeys(ele,bufk);
    int nbdofs=bufk.size();
    int nbcomp=comp.size();
    int curpos=keys.size();
    keys.reserve(curpos+nbcomp*nbdofs);
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
};

class DgC0FullDgLagrangeFunctionSpace : public DgC0LagrangeFunctionSpace, public DgC0FullDgFunctionSpaceBase{
 public :
  DgC0FullDgLagrangeFunctionSpace(int id) : DgC0LagrangeFunctionSpace(id){};
  ~DgC0FullDgLagrangeFunctionSpace(){};

  virtual void getKeys(MElement *ele, std::vector<Dof> &keys){
    int nk=ScalarFS->getNumKeys(ele); // return the number of vertices
    std::vector<Dof> bufk;
    bufk.reserve(nk);
    ScalarFS->getKeys(ele,bufk);
    int nbdofs=bufk.size();
    int nbcomp=comp.size();
    int curpos=keys.size();
    keys.reserve(curpos+nbcomp*nbdofs);
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

  virtual void getKeys(MElement *ele, std::vector<Dof> &keys, std::vector<int> *numvertex)
  {
    int nk=ScalarFS->getNumKeys(ele); // return the number of vertices
    std::vector<Dof> bufk;
    bufk.reserve(nk);
    ScalarFS->getKeys(ele,bufk);
    int nbdofs=bufk.size();
    int nbcomp=comp.size();
    int curpos=keys.size();
    keys.reserve(curpos+nbcomp*nbdofs);
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
};

#endif // DGC0PLATEFUNCTIONSPACE_H_
