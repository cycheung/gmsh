//
// C++ Interface: terms
//
// Description: Function to Assemble terms for plate element (bulk, interface and interface boundary
//              TODO Template this function with Solver/solverAlgorithms
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef DGC0PLATEALGORITHMS_H_
#define DGC0PLATEALGORITHMS_H_
#include "dofManager.h"
#include "DgC0PlateTerms.h"
#include "quadratureRules.h"
#include "MVertex.h"
#include "MInterfaceElement.h"
#include "DgC0PlateFunctionSpace.h"

// One !=ce Dynamic cast on function space to have the good getKeysFunction oninterface
// and iterator on interfaceElement
template<class Iterator,class Assembler> void MyAssembleOnInterface(BilinearTermBase &term,FunctionSpaceBase &space,Iterator itbegin,Iterator itend,QuadratureBase &integrator,
                                                         Assembler &assembler) // symmetric
{
  fullMatrix<typename Assembler::dataMat> localMatrix;
  std::vector<Dof> R;
  DgC0FunctionSpaceBase *dgspace = dynamic_cast<DgC0FunctionSpaceBase*>(&space);
  for (Iterator it = itbegin;it!=itend; ++it)
  {
    MInterfaceElement *e = *it;
    R.clear();
    IntPt *GP;
    int npts=integrator.getIntPoints(e,&GP);
    term.get(e,npts,GP,localMatrix);
    dgspace->getKeys(e,R);
    assembler.assemble(R, localMatrix);
  }
}
template<class Iterator,class Assembler> void MyAssembleOnInterface(LinearTermBase &term,FunctionSpaceBase &space,Iterator itbegin,
                                                         Iterator itend,QuadratureBase &integrator,Assembler &assembler) // symmetric
{
  fullVector<typename Assembler::dataMat> localVector;
  std::vector<Dof> R;
  DgC0FunctionSpaceBase *dgspace = dynamic_cast<DgC0FunctionSpaceBase*>(&space);
  for (Iterator it = itbegin;it!=itend; ++it)
  {
    MInterfaceElement *e = *it;
    R.clear();
    IntPt *GP;
    int npts=integrator.getIntPoints(e,&GP);
    term.get(e,npts,GP,localVector);
    dgspace->getKeys(e,R);
    assembler.assemble(R, localVector);
  }
}


template<class Iterator,class Assembler> void Assemble(LinearTermBase &term,FunctionSpaceBase &space,Iterator itbegin,
                                                       Iterator itend,QuadratureBase &integrator,Assembler &assembler,
                                                       std::vector<MInterfaceElement*> vinter)
{
  fullVector<typename Assembler::dataMat> localVector;
  std::vector<Dof> R;
  for (Iterator it = itbegin;it!=itend; ++it)
  {
    MElement *e = *it;
    MInterfaceElement *ielem=NULL;
    int ind,kind;
    // Find the interfaceElement associated to the element where the boundary condition is applied
    // The vertex of *e chosen for test depends on the BC is on a vertex or edge If condition is applied on a line
    if(e->getNumVertices()==1) ind=0; //test on neumann.onWhat ??
    else ind=2;
      for(std::vector<MInterfaceElement*>::iterator it1 = vinter.begin(); it1!=vinter.end();++it1){
        int nn=(*it1)->getNumVertices();
        for(int k=ind;k<nn;k++)
          if (e->getVertex(ind)==(*it1)->getVertex(k)) {ielem=*it1;kind=k;break;} // kind used only for BC on vertex
       }
    if(ielem==NULL) printf("Warning : impossible to applied the Neumann Boundary Condition on element %d",e->getNum());
    // node where the boundary conditions has to be applied
    int nvi=ielem->getNumVertices();
    std::vector<int> vernum;
    std::vector<int> te;
    vernum.reserve(nvi);
    for(int i=0;i<nvi;i++) vernum.push_back(0);
      ielem->getLocalVertexNum(0,vernum);
    std::vector<int> *pver = &vernum;
    if(ind==0){
      te.push_back(vernum[kind]);
      pver = &te;
    }
    R.clear();
    IntPt *GP;
    int npts=integrator.getIntPoints(e,&GP);
    term.get(e,npts,GP,localVector);
    DgC0FullDgFunctionSpaceBase* dgspace = dynamic_cast<DgC0FullDgFunctionSpaceBase*>(&space);
    dgspace->getKeys(ielem->getElem(0),R,pver);
    assembler.assemble(R, localVector);
  }
}

// Dirichlet
// TODO regroup this fuinction with FixNodalDofs
template<class Assembler> void FixNodalDofs(DgC0FullDgFunctionSpaceBase &space,MElement *e,Assembler &assembler,
                                            simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,
                                            std::vector<MInterfaceElement*> &inter, std::vector<MInterfaceElement*> &internalInterface)
{
  bool InternalEdge=false;
  std::vector<MVertex*> tabV;
  int nv=e->getNumVertices();
  for (int i=0;i<nv;++i) tabV.push_back(e->getVertex(i));
  MInterfaceElement *ielem=NULL;
  // Find the interface element linked with e
  int ind, posfind;
  if(nv==1) ind =0; // BC to applied on a vertex --> look at extremities too
  else ind = 2;
  for(std::vector<MInterfaceElement*>::iterator it=inter.begin();it!=inter.end();++it){
    int nn = (*it)->getNumVertices();
    for(int i=ind;i<nn;i++)  // pass node on extremities
      if(tabV[ind]==(*it)->getVertex(i)) {ielem=*it; posfind = i;}
  }
  // If not find on external edge. Look in internal edge
  if(ielem==NULL){
    for(std::vector<MInterfaceElement*>::iterator it=internalInterface.begin();it!=internalInterface.end();++it){
      int nn = (*it)->getNumVertices();
      for(int i=ind;i<nn;i++)  // pass node on extremities
        if(tabV[ind]==(*it)->getVertex(i)) {ielem=*it; posfind = i; InternalEdge=true;}
    }
  }
  std::vector<Dof> R;
  if(ielem==NULL) printf("Impossible to fix dof for element %d\n",e->getNum());
  else{
    int nvi=ielem->getNumVertices();
    std::vector<int> vernum;
    vernum.resize(nvi);
    //for(int i=0;i<nvi;i++) vernum.push_back(0);
    ielem->getLocalVertexNum(0,vernum);
    if(nv == 1){ // if boundary condition is applied on a vertex
      int temp = vernum[posfind];
      vernum.clear();
      vernum.push_back(temp);
    }
    std::vector<int> *pver = &vernum;
//    DgC0FunctionSpaceBase* dgspace = dynamic_cast<DgC0FunctionSpaceBase*>(&space);
    space.getKeys(ielem->getElem(0),R,pver);
    if(InternalEdge){
      ielem->getLocalVertexNum(1,vernum);
      if(nv == 1){ // if boundary condition is applied on a vertex
        int temp = vernum[posfind];
        vernum.clear();
        vernum.push_back(temp);
      }
      space.getKeys(ielem->getElem(1),R,pver);
    }
    for (std::vector<Dof>::iterator itd=R.begin();itd!=R.end();++itd)
    {
      Dof key=*itd;
      if (filter(key))
      {
        for (int i=0;i<nv;++i)
        {
          if ((ielem->getElem(0)->getNum()==key.getEntity()) or (ielem->getElem(1)->getNum() == key.getEntity()))
          {
            //printf("Dof %d %d fixed\n",key.getEntity(),key.getType());
            assembler.fixDof(key, fct(tabV[i]->x(),tabV[i]->y(),tabV[i]->z()));
            break;
          }
        }
      }
    }
  }
}

template<class Iterator,class Assembler> void FixNodalDofs(DgC0FullDgFunctionSpaceBase &space,Iterator itbegin,Iterator itend,Assembler &assembler,
                                                           simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,
                                                           std::vector<MInterfaceElement*> &inter,std::vector<MInterfaceElement*> &vinternalinter)
{
  for (Iterator it=itbegin;it!=itend;++it)
  {
    FixNodalDofs(space,*it,assembler,fct,filter,inter,vinternalinter);
  }
}

// Allow to know the total number of dof
template<class Iterator,class Assembler> void NumberDofs(FunctionSpaceBase &space,Iterator itbegin,Iterator itend,Assembler &assembler,bool fullDg)
{
 for (Iterator it=itbegin;it!=itend;++it)
  {
    MElement *e=*it;
    std::vector<Dof> R;
    DgC0FunctionSpace<SVector3> *dgspace = dynamic_cast<DgC0FunctionSpace<SVector3>*>(&space);
    dgspace->getKeys(e,R);
    int nbdofs=R.size();
    for (int i=0;i<nbdofs;++i) assembler.numberDof(R[i]);
  }
}

template<class Assembler> void FixNodalDofs(FunctionSpaceBase &space,MElement *e,Assembler &assembler,simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,bool fullDg)
{
  std::vector<MVertex*> tabV;
  int nv=e->getNumVertices();
  std::vector<Dof> R;
  DgC0FunctionSpace<SVector3> *dgspace = dynamic_cast<DgC0FunctionSpace<SVector3>*>(&space);
  dgspace->getKeys(e,R);
  tabV.reserve(nv);
  for (int i=0;i<nv;++i) tabV.push_back(e->getVertex(i));

  if(!fullDg){
    for (std::vector<Dof>::iterator itd=R.begin();itd!=R.end();++itd)
    {
      Dof key=*itd;
      if (filter(key))
      {
        for (int i=0;i<nv;++i)
        {
          if (tabV[i]->getNum()==key.getEntity())
          {
            //printf("Fix dof number %d comp %d\n",key.getEntity(),key.getType());
            assembler.fixDof(key, fct(tabV[i]->x(),tabV[i]->y(),tabV[i]->z()));
            break;
          }
        }
      }
    }
  }
  else{
    for (std::vector<Dof>::iterator itd=R.begin();itd!=R.end();++itd)
    {
      Dof key=*itd;
      if (filter(key))
      {
        for (int i=0;i<nv;++i)
        {
          if (e->getNum()==key.getEntity())
          {
            //printf("Fix dof number %d comp %d\n",key.getEntity(),key.getType());
            assembler.fixDof(key, fct(tabV[i]->x(),tabV[i]->y(),tabV[i]->z()));
            break;
          }
        }
      }
    }
  }
}

template<class Iterator,class Assembler> void FixNodalDofs(FunctionSpaceBase &space,Iterator itbegin,Iterator itend,Assembler &assembler,simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,bool fullDg)
{
  for (Iterator it=itbegin;it!=itend;++it)
  {
    FixNodalDofs(space,*it,assembler,fct,filter,fullDg);
  }
}

template<class Iterator,class Assembler> void Assemble(LinearTermBase &term,FunctionSpaceBase &space,Iterator itbegin,Iterator itend,QuadratureBase &integrator,
                                                       Assembler &assembler,bool fullDg)
{
  fullVector<typename Assembler::dataMat> localVector;
  std::vector<Dof> R;
  for (Iterator it = itbegin;it!=itend; ++it)
  {
    MElement *e = *it;
    R.clear();
    IntPt *GP;
    int npts=integrator.getIntPoints(e,&GP);
    term.get(e,npts,GP,localVector);
    DgC0FunctionSpace<SVector3> *dgspace = dynamic_cast<DgC0FunctionSpace<SVector3>*>(&space);
    dgspace->getKeys(e,R);
    assembler.assemble(R, localVector);
  }
}

class DgC0PlateFilterDofComponent :public FilterDof
{
  int comp;
 public :
  DgC0PlateFilterDofComponent(int comp_): comp(comp_) {}
  virtual bool operator()(Dof key)
  {
    int type=key.getType();
    int icomp,iphys,inode;
    DgC0PlateDof::getThreeIntsFromType(type, icomp, iphys,inode);
    if (icomp==comp) return true;
    return false;
  }
};

#endif // DGC0PLATEALGORITHMS_H_
