//
// C++ Interface: terms
//
// Description: Solver field for C0Dg plate
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef _DGC0PLATESOLVERFIELD_H_
#define _DGC0PLATESOLVERFIELD_H_

#include <vector>
#include <iostream>
#include "MElement.h"
#include "dofManager.h"
#include "DgC0PlateFunctionSpace.h"

template<class T>
class DgC0SolverField //: public DgC0FunctionSpace<T> // being able to use it instead of a real function space is interesting (nbkeys=1, explicit keys/dofs undefined (or could be defined element-wise )
{
 public:
  typedef typename TensorialTraits<T>::ValType ValType;
  typedef typename TensorialTraits<T>::GradType GradType;
  typedef typename TensorialTraits<T>::HessType HessType;
 private:
  dofManager<double> *dm;
  DgC0FunctionSpace<T> *fs;
 public:
  DgC0SolverField(dofManager<double> *dm_,  DgC0FunctionSpace<T> *fs_) : dm(dm_),fs(fs_) {}
  virtual int getNumKeys(MVertex *ver) { return 1;}
  virtual int getNumKeys(MElement *ele) { return 1;}
 private:
  virtual void getKeys(MElement *ele, std::vector<Dof> &keys) { Msg::Error("getKeys for SolverField should'nt be called");}
  virtual void getKeys(MVertex *ver, std::vector<Dof> &keys) {Msg::Error("getKeys for SolverField should'nt be called");}
 public:

/*  virtual void f(MElement *ele, double u, double v, double w, ValType &val)
  {
    int field,comp;
    std::vector<Dof> D;
    std::vector<double> SFVals;
    std::vector<double> DMVals;
    fs->getKeys(ele,D,false);
    dm->getDofValue(D,DMVals);
    fs->f(ele,u,v,w,SFVals);
    val=ValType();
    for (int i=0;i<D.size();i++){
        D[i].getTwoIntsFromType(D[i].getType(), comp,field);
        val(comp)+=SFVals[0]*DMVals[i];
    }
  }*/

  virtual void getVertexDisplacement(MElement *ele,ValType &val,bool FullDg, int vernum=0) // vernum used for fullDg
  {
    int ent;
    if(!FullDg) ent=ele->getVertex(0)->getNum();
    else ent = ele->getNum();
    double temp;
    for(int j=0;j<3;j++){ // Pass number of component and loop on it
      Dof D(ent,DgC0PlateDof::createTypeWithThreeInts(j,1000,vernum)); // create a dof number getId function ??
      //temp = dm->getDofValue(D);
      dm->getDofValue(D,temp);
      val(j) = temp;
    }
  }

/*  virtual void f(MElement *ele, double u, double v, double w, std::vector<ValType> &vals)
  {
    ValType val;
    f(ele,u,v,w,val);
    vals.push_back(val);
  }
*/
/*  virtual void gradf(MElement *ele, double u, double v, double w,GradType &grad)
  {
    int field,comp;
    std::vector<Dof> D;
    std::vector<SVector3> SFGrads;
    std::vector<double> DMVals;
    fs->getKeys(ele,D);
    dm->getDofValue(D,DMVals);
    fs->gradf(ele,u,v,w,SFGrads);
    grad=GradType();
    int n = ele->getNumVertices();
    for (int i=0;i<D.size();i++){
      D[i].getTwoIntsFromType(D[i].getType(), comp, field);
      for(int j=0;j<3;j++)
        for(int k=0;k<n;k++)
          grad(j,comp)+= SFGrads[k][j] * DMVals[i];
    }
  }


  virtual void gradf(MElement *ele, double u, double v, double w,std::vector<GradType> &grads)
  {
    GradType grad;
    gradf(ele,u,v,w,grad);
    grads.push_back(grad);
  }
    virtual void hessfuvw(MElement *ele, double u, double v, double w,std::vector<HessType> &hess)
  {
    //HessType hes;
    fs->hessfuvw(ele,u,v,w,hess);
  }*/
};

#endif //_DGC0PLATESOLVERFIELD_H_

