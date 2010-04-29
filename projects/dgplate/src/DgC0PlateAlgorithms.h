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

// To not pollute Solver/solverAlgorithms but can be put into this file
template<class Iterator,class Assembler> void AssembleBulk(DgC0BilinearTermBase &term,DgC0FunctionSpaceBase &space,Iterator itbegin,Iterator itend,QuadratureBase &integrator,Assembler &assembler, bool FullDg) // symmetric
{
  fullMatrix<typename Assembler::dataMat> localMatrix;
  std::vector<Dof> R;
  for (Iterator it = itbegin;it!=itend; ++it)
  {
    MElement *e = *it;
    R.clear();
    IntPt *GP;
    int npts=integrator.getIntPoints(e,&GP);
    term.get(e,npts,GP,localMatrix);
    int nv=e->getNumVertices();
    std::vector<MVertex*> ver;
    ver.reserve(nv);
    for(int i=0;i<nv;i++) ver.push_back(e->getVertex(i));
    space.getKeys(e,R,FullDg);
    assembler.assemble(R, localMatrix);
  }
}

template<class Assembler> void AssembleInterface(DgC0BilinearTermBase &iterm,DgC0FunctionSpaceBase &space, groupOfElements* ve, std::vector<MInterfaceElement*> &vie, QuadratureBase &integrator,Assembler &assembler, bool FullDg){
  // creation of a matrix this matrix is resize and put = 0 in the function iterm.get()
  fullMatrix<double> localMatrix;
  // vector with dof the dof are append in space.getKeys()
  std::vector<Dof> R;
  // Loop on interface element
  for (int i = 0; i < vie.size(); i++)
  {
    MInterfaceElement *ie = vie[i];
    // get pointers on the elements linked to the interface element i
    MElement **elem = ie->getElem();
    // Clear the dof vector (has to be clear because getKeys is an append method)
    R.clear();
    // Compute the integration point
    IntPt *GP;
    int npts=integrator.getIntPoints(ie,&GP);
    // compute elementatry matrix
    iterm.get(ie,npts,GP,localMatrix);
    // dof
    if(ie->getElem(1) == ie->getElem(0)){ // Element is on the boundary
      space.getKeys(ie->getElem(0),R,FullDg);
    }
    else{
      space.getKeys(ie->getElem(0),R,FullDg); space.getKeys(ie->getElem(1),R,FullDg);
    }
  assembler.assemble(R, localMatrix);
  }
}

// Dirichlet
// TODO regroup this fuinction with FixNodalDofs
template<class Assembler> void FixNodalDofs(DgC0FunctionSpaceBase &space,MElement *e,Assembler &assembler,simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,std::vector<MInterfaceElement*> &inter)
{
  std::vector<MVertex*> tabV;
  int nv=e->getNumVertices();
  tabV.reserve(nv);
  for (int i=0;i<nv;++i) tabV.push_back(e->getVertex(i));
  MInterfaceElement *ielem=NULL;
  // Find the interface element linked with e
  for(std::vector<MInterfaceElement*>::iterator it=inter.begin();it!=inter.end();++it){
    int nn = (*it)->getNumVertices();
    for(int i=2;i<nn;i++)  // pass node on extremities
      if(tabV[2]==(*it)->getVertex(i)) ielem=*it;
  }
  std::vector<Dof> R;
  if(ielem==NULL) printf("Impossible to fix dof for element %d\n",e->getNum());
  else{
    int nvi=ielem->getNumVertices();
    std::vector<int> vernum;
    vernum.reserve(nvi);
    for(int i=0;i<nvi;i++) vernum.push_back(0);
    ielem->getLocalVertexNum(0,vernum);
    std::vector<int> *pver = &vernum;
    space.getKeys(ielem,R,true,pver);
    for (std::vector<Dof>::iterator itd=R.begin();itd!=R.end();++itd)
    {
      Dof key=*itd;
      if (filter(key))
      {
        for (int i=0;i<nv;++i)
        {
          if (ielem->getElem(0)->getNum()==key.getEntity()) // TODO move 1e5 in dofManager
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

template<class Iterator,class Assembler> void FixNodalDofs(DgC0FunctionSpaceBase &space,Iterator itbegin,Iterator itend,Assembler &assembler,simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,std::vector<MInterfaceElement*> &inter)
{
  for (Iterator it=itbegin;it!=itend;++it)
  {
    FixNodalDofs(space,*it,assembler,fct,filter,inter);
  }
}

template<class Iterator,class Assembler> void NumberDofs(DgC0FunctionSpaceBase &space,Iterator itbegin,Iterator itend,Assembler &assembler,bool fullDg)
{
 for (Iterator it=itbegin;it!=itend;++it)
  {
    MElement *e=*it;
    std::vector<Dof> R;
    space.getKeys(e,R,fullDg);
    int nbdofs=R.size();
    for (int i=0;i<nbdofs;++i) assembler.numberDof(R[i]);
  }
}

template<class Iterator,class Assembler> void Assemble(DgC0LinearTermBase &term,DgC0FunctionSpaceBase &space,Iterator itbegin,Iterator itend,QuadratureBase &integrator,Assembler &assembler, std::vector<MInterfaceElement*> vinter)
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
    vernum.reserve(nvi);
    for(int i=0;i<nvi;i++) vernum.push_back(0);
      ielem->getLocalVertexNum(0,vernum);
    std::vector<int> *pver = &vernum;
    if(ind==0){
      std::vector<int> te;
      te.push_back(vernum[kind]);
      pver = &te;
    }
    R.clear();
    IntPt *GP;
    int npts=integrator.getIntPoints(e,&GP);
    term.get(e,npts,GP,localVector);
    space.getKeys(ielem->getElem(0),R,true,pver);
    assembler.assemble(R, localVector);
  }
}

template<class Assembler> void FixNodalDofs(DgC0FunctionSpaceBase &space,MElement *e,Assembler &assembler,simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,bool fullDg)
{
  std::vector<MVertex*> tabV;
  int nv=e->getNumVertices();
  std::vector<Dof> R;
  space.getKeys(e,R,fullDg);
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

template<class Iterator,class Assembler> void FixNodalDofs(DgC0FunctionSpaceBase &space,Iterator itbegin,Iterator itend,Assembler &assembler,simpleFunction<typename Assembler::dataVec> &fct,FilterDof &filter,bool fullDg)
{
  for (Iterator it=itbegin;it!=itend;++it)
  {
    FixNodalDofs(space,*it,assembler,fct,filter,fullDg);
  }
}

template<class Iterator,class Assembler> void Assemble(DgC0LinearTermBase &term,DgC0FunctionSpaceBase &space,Iterator itbegin,Iterator itend,QuadratureBase &integrator,Assembler &assembler,bool fullDg)
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
    space.getKeys(e,R,fullDg);
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
