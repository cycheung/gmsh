//
// C++ Interface: terms
//
// Description: dofManager for non linear system (rewritte of assemble(fullMatrix) no contribution for fixed Dof
//              The contribution is taken into account in the assemble(fullvector) of RightHandSide
//              More add functions to archiving force can be placed in dofManager but I put it here to not pollute gmsh code
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef DOFMANAGERNONLINEARSYSTEM_H_
#define DOFMANAGERNONLINEARSYSTEM_H_

#include "dofManager.h"
template<class T>
class dofManagerNLS : public dofManager<T>{
 private :
  // map to retains the force of fixed dof (needed for archiving) the other RHS value are accessible in RHS
  std::map<Dof, typename dofManager<T>::dataVec> RHSfixed;

 public :

 dofManagerNLS(linearSystem< typename dofManager<T>::dataMat > *l,
               std::vector<Dof> varchdof, typename dofManager<T>::dataVec value) : dofManager<T>(l){
    for(int i=0;i<varchdof.size();i++)
      RHSfixed[varchdof[i]] = value;

}
 dofManagerNLS(linearSystem< typename dofManager<T>::dataMat > *l1, linearSystem<typename dofManager<T>::dataMat> *l2,
               std::vector<Dof> &varchdof, typename dofManager<T>::dataVec &value) : dofManager<T>(l1,l2){
    for(int i=0;i<varchdof.size();i++)
      RHSfixed[varchdof[i]] = value;
}

 // this function is replaced to avoid a double contribution in right hand side for fixed dof
  inline void assemble(std::vector<Dof> &R, const fullMatrix<typename dofManager<T>::dataMat> &m)
  {
    if (!this->_current->isAllocated()) this->_current->allocate(this->unknown.size());
    std::vector<int> NR(R.size());
    for (unsigned int i = 0; i < R.size(); i++)
    {
      std::map<Dof, int>::iterator itR = this->unknown.find(R[i]);
      if (itR != this->unknown.end()) NR[i] = itR->second;
      else NR[i] = -1;
    }
    for (unsigned int i = 0; i < R.size(); i++)
    {
      if (NR[i] != -1)
      {
        for (unsigned int j = 0; j < R.size(); j++)
        {
          if (NR[j] != -1)
          {
            this->_current->addToMatrix(NR[i], NR[j], m(i, j));
          }
          else
          {
            typename std::map<Dof,  typename dofManager<T>::dataVec>::iterator itFixed = this->fixed.find(R[j]);
            if (itFixed == this->fixed.end())
              assembleLinConst(R[i],R[j],m(i,j));
          }
        }
      }
      else
      {
        for (unsigned int j = 0; j < R.size(); j++){
        assembleLinConst(R[i],R[j],m(i,j));
        }
      }
    }
  }
  // initial value is needed because dataVect is unknown FIX IT HOW ??
  void archivingRHS(Dof key, typename dofManager<T>::dataVec &value){
    std::map<Dof, int>::iterator itR = this->unknown.find(key);
      if (itR == this->unknown.end()) RHSfixed[key]=value; // otherwise the value is accessible in RHS vector
  }
  // function assemble(vect) is rewritten to archiving rhs if needed
  virtual inline void assemble(std::vector<Dof> &R, const fullVector<typename dofManager<T>::dataMat> &m) // for linear forms
  {
    if (!this->_current->isAllocated()) this->_current->allocate(this->unknown.size());
    std::vector<int> NR(R.size());
    // Archiving
    for (unsigned int i = 0; i < R.size(); i++)
    {
      typename std::map<Dof,  typename dofManager<T>::dataVec>::iterator itR = RHSfixed.find(R[i]);
      if (itR != RHSfixed.end()) itR->second += m(i); // add of contribution set to zero when the value is by a function FIX it
      else NR[i] = -1;
    }

    for (unsigned int i = 0; i < R.size(); i++)
    {
      std::map<Dof, int>::iterator itR = this->unknown.find(R[i]);
      if (itR != this->unknown.end()) NR[i] = itR->second;
      else NR[i] = -1;
    }
    for (unsigned int i = 0; i < R.size(); i++)
    {
      if (NR[i] != -1)
      {
        this->_current->addToRightHandSide(NR[i], m(i));
      }
      else
      {
        typename std::map<Dof,DofAffineConstraint<typename dofManager<T>::dataVec> >::iterator itConstraint;
        itConstraint = this->constraints.find(R[i]);
        if (itConstraint != this->constraints.end())
        {
          for (unsigned i=0;i<(itConstraint->second).linear.size();i++)
          {
                  std::map<Dof, int>::iterator itC = this->unknown.find((itConstraint->second).linear[i].first); // lin dep in unknown ?!
                  this->_current->addToRightHandSide(itC->second, m(i)*(itConstraint->second).linear[i].second);
          }
        }
      }
    }
  }
  // zero for RHSfixed
  void clearRHSfixed(){
    for(typename std::map<Dof,typename dofManager<T>::dataVec>::iterator itR = RHSfixed.begin();itR != RHSfixed.end(); ++itR)
      itR->second = 0.;
  }

  void getForces(std::map<int,std::vector<Dof> > &aef, std::map<int,typename dofManager<T>::dataVec> &data){
    for(std::map<int,std::vector<Dof> >::iterator it = aef.begin(); it != aef.end(); ++it){
      double Fint =0.;
      for(int i=0;i<it->second.size();i++)
        Fint -= RHSfixed[it->second[i]]; // as internal forces are computed with a minus sign in computeRightHandSide;
    data[it->first] = Fint;
    }
  }

};

#endif // DOFMANAGERNONLINEARSYSTEM_H_
