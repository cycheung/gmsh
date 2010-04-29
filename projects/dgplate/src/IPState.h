//
// C++ Interface: terms
//
// Description: Class to store internal variables at gauss point
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef IPSTATE_H_
#define IPSTATE_H_

#include "GModelWithInterface.h"
#include "DgC0PlateSolver.h"
#include "MInterfaceElement.h"
#include "MElement.h"
#include "LocalBasis.h"
#include "materialLaw.h"
#include "reduction.h"
#include "DgC0PlateSolverField.h"
// enum to access to component of stress and deformation GLOBAL ??
struct component{
 enum enumcomp{xx,yy,zz,xy,xz,yz};
};

// class to number to IP
class IPnum{
  protected :
    long int _entity;
    int _type;
  public :
    IPnum(long int ent, int typ) : _entity(ent),_type(typ){};
    inline static int createTypeWithTwoInts(const int i1, const int i2){
      return 1000*i1+i2;
    }
    void getTypeFromTwoInt(const int t, int &i1, int &i2){
      i1 = t/1000;
      i2 = t%1000;
    }
    std::pair<long int, int> get() const {
      return std::pair<long int,int>(_entity,_type);
    }
    long int getEntity() const{return _entity;}
    int getType() const{return _type;}
};

// class with the variables of IP (stress, deformation and localBasis)
class IPVariable
{
  public :
    IPVariable(){}
    ~IPVariable(){}
    // copie constructor
    IPVariable(const IPVariable &source){
    }
    IPVariable &operator = (const IPVariable &source){
    }
    virtual void bidon(){};

};
class IPVariablePlate : public IPVariable{
  protected :
    double sigma[6];
    double epsilon[6];
    LocalBasis lb;
  public :
    void computeStressAndDeformation(materialLaw *mtl,const int nbFF, const int nbdof, const fullMatrix<double> &disp, const std::vector<TensorialTraits<double>::GradType> &Grads, const std::vector<TensorialTraits<double>::HessType> &Hess){
      // Deformation (small deformation plate for now)
      epsilon[0] = epsilongd(0,0,&lb,Grads,0,disp) + rhogd(0,0,&lb,Hess,0,disp); // TODO rewrite epsilongd and rhogd more efficiently
      epsilon[1] = epsilongd(1,1,&lb,Grads,0,disp) + rhogd(1,1,&lb,Hess,0,disp);
      epsilon[3] = epsilongd(0,1,&lb,Grads,0,disp) + rhogd(0,1,&lb,Hess,0,disp);
      epsilon[2] = epsilon[4] = epsilon[5] = 0.;
      // stress thanks to material law
      mtl->stress(epsilon,sigma);
      // multiplication by jacobian ?? (in stress function ??)
      double invJ=1./lb.getJacobian();
      for(int i=0;i<6;i++) sigma[i]*=invJ;
    }
    // On interface the local basis of Element is used to compute and it's given as argument
    void computeStressAndDeformation(materialLaw *mtl,const LocalBasis *lbe,const int nbFF, const int nbdof, const fullMatrix<double> &disp, const std::vector<TensorialTraits<double>::GradType> &Grads, const std::vector<TensorialTraits<double>::HessType> &Hess){
      // Deformation (small deformation plate for now)
      epsilon[0] = epsilongd(0,0,lbe,Grads,0,disp) + rhogd(0,0,lbe,Hess,0,disp); // TODO rewrite epsilongd and rhogd more efficiently
      epsilon[1] = epsilongd(1,1,lbe,Grads,0,disp) + rhogd(1,1,lbe,Hess,0,disp);
      epsilon[3] = epsilongd(0,1,lbe,Grads,0,disp) + rhogd(0,1,lbe,Hess,0,disp);
      epsilon[2] = epsilon[4] = epsilon[5] = 0.;
      // stress thanks to material law
      mtl->stress(epsilon,sigma);
      // multiplication by jacobian ?? (in stress function ??)
      double invJ=1./lbe->getJacobian();
      for(int i=0;i<6;i++) sigma[i]*=invJ;
    }
    // As template constructor doesn't take any argument (A set function allow to compute initial data)
    IPVariablePlate() : IPVariable() {for(int i=0;i<6;i++){sigma[i]=0.;epsilon[i]=0.;}}
    ~IPVariablePlate(){}
    // copie constructor
    IPVariablePlate(const IPVariablePlate &source):IPVariable(source)
    {
      for(int i=0;i<6;i++){
        sigma[i]=source.getSigma(i);
        epsilon[i]=source.getEpsilon(i);
      }
      lb=source.lb;
    }
    IPVariablePlate &operator = (const IPVariablePlate &source){
      IPVariable::operator = (source);
      for(int i=0;i<6;i++){
        sigma[i]=source.getSigma(i);
        epsilon[i]=source.getEpsilon(i);
      }
      lb=source.lb;
      return *this;
    }
    virtual double getSigma(const int i) const{return sigma[i];}
    virtual double getEpsilon(const int i) const{return epsilon[i];}
    virtual void setLocalBasis(MElement *e, const std::vector<TensorialTraits<double>::GradType> &Grads){
      lb.set(e,Grads);
    }
    virtual void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads, const SVector3 &t0p, const SVector3 &t0m)
    {
      lb.set(ie,Grads,t0p,t0m);
    }
    virtual void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads, const SVector3 &t0m){
      lb.set(ie,Grads,t0m);
    }
};

class IPState{
  protected :
    IPVariable *_initial;     // initial state t=0
    IPVariable *_previous;    // previous step t=n-1
    IPVariable *_current;     // current step
    // pointer on a bool value to choice what vector is current and what vector is previous
    bool *_st;
    // pointer to the elasticField of IP
    DGelasticField* _ef;
  public :
    enum whichState{initial, previous, current}; //protected enum ??
    IPState(bool *st,DGelasticField* ef) : _st(st), _ef(ef)
    {
      // should be done by material trype
      _initial = new IPVariablePlate();
      _previous = new IPVariablePlate();
      _current = new IPVariablePlate();
    }
    ~IPState(){
        if(_initial) delete _initial;
        if(_current) delete _current;
        if(_previous) delete _previous;

    }
    IPState(const IPState &source){
        // should be done by material trype
        _initial = new IPVariablePlate();
        _previous = new IPVariablePlate();
        _current = new IPVariablePlate();

      _initial->operator = (*(source._initial));
      _previous = source._previous;
      _current = source._current;
      _st = source._st;
      _ef=source._ef;
    }
    IPState & operator = (const IPState &source){
      _initial->operator = (*(source._initial));
      _initial = source._initial;
      _previous = source._previous;
      _current = source._current;
      _st = source._st;
      _ef=source._ef;
      return *this;
    }
    IPVariable* getState(const whichState wst) const {
      switch(wst){
        case initial :
          return _initial;
          break;
        case previous :
          if(*_st) return _previous; else return _current;
          break;
        case current :
          if(*_st) return _current; else return _previous;
          break;
        default : Msg::Error("Impossible to select the desired state for internal variable \n");
      }
    }


 /*   void computeStressAndDeformation(whichState ws,materialLaw *mtl,const int nbFF, const int nbdof, const fullMatrix<double> &disp, const std::vector<TensorialTraits<double>::GradType> &Grads, const std::vector<TensorialTraits<double>::HessType> &Hess){
      switch(ws){
        case initial :
          _initial->computeStressAndDeformation(mtl,nbFF,nbdof,disp,Grads,Hess);
          break;
        case previous :
          _previous->computeStressAndDeformation(mtl,nbFF,nbdof,disp,Grads,Hess);
          break;
        case current :
          _current->computeStressAndDeformation(mtl,nbFF,nbdof,disp,Grads,Hess);
          break;
        default : Msg::Error("Impossible to select the desired state for internal variable \n");
      }
    }*/
};

// Class to access to the IPState of all gauss point
class AllIPState{
  protected :
    std::vector < IPState > allIS;
    typedef std::map<std::pair<long int,int>, IPState*> ipstateContainer;
    ipstateContainer _mapall;
    // flag to switch previous and current (change the pointer in place of copy all variables)
    bool state;
    void createAndStoreIP(IPnum *key,DGelasticField *ef){
      //IPState ipstatmp(&state,ef);
      //allIS.push_back(ipstatmp);
      //std::pair<long int, int> tp(key->getEntity(),key->getType());
      //_mapall.insert(std::pair<std::pair<long int,int>,IPState*>(tp,&allIS.back()));
      IPState*ips = new IPState(&state,ef);
      //allIS.push_back(ipstatmp);
      std::pair<long int, int> tp(key->getEntity(),key->getType());
      _mapall.insert(std::pair<std::pair<long int,int>,IPState*>(tp,ips));
    }

  public :
    // iterator ATTENTION cannot definied with the template T so must be define for each IPvariable
    std::map<std::pair<long int,int>, IPState*>::const_iterator begin() const {return _mapall.begin();}
    std::map<std::pair<long int,int>, IPState*>::const_iterator end() const {return _mapall.end();}

    AllIPState(GModelWithInterface *pModel, std::vector<DGelasticField> elasf, QuadratureBase &integrator_on_bulk, QuadratureBase &integrator_on_inter){
      state = true; // at creation of object state is true
      IntPt *GP; // needed to know the number of gauss point
      for(int i=0;i<elasf.size();i++){
        if(!elasf[i].getFormulation()){ //cg/dg formulation
          // The object is initialized thanks to interfaceElement. Other Possibilities ?? if CG can be done with element but not for cg/dg and full dg
          // loop
          for(std::vector<MInterfaceElement*>::iterator it=elasf[i].gi.begin(); it!=elasf[i].gi.end();++it){
            MInterfaceElement *ie = *it;
            // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
            MElement *em = ie->getElem(0);
            int edge = ie->getEdgeNumber(0);
            int npts_inter=integrator_on_inter.getIntPoints(ie,&GP);
            for(int j=0;j<npts_inter;j++){
              IPnum key=IPnum(em->getNum(),IPnum::createTypeWithTwoInts(edge,j));
              createAndStoreIP(&key,&elasf[i]);
            }
          }
        }
        else{ //full dg formulation
          // loop
          for(std::vector<MInterfaceElement*>::iterator it=elasf[i].gi.begin(); it!=elasf[i].gi.end();++it){
            MInterfaceElement *ie = *it;
            // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
            MElement *em = ie->getElem(0);
            MElement *ep = ie->getElem(1);
            int edgem = ie->getEdgeNumber(0);
            int edgep = ie->getEdgeNumber(1);
            int npts_inter=integrator_on_inter.getIntPoints(ie,&GP);
            for(int j=0;j<npts_inter;j++){
              IPnum key=IPnum(em->getNum(),IPnum::createTypeWithTwoInts(edgem,j));
              createAndStoreIP(&key,&elasf[i]);
              IPnum keyp=IPnum(ep->getNum(),IPnum::createTypeWithTwoInts(edgep,j));
              createAndStoreIP(&keyp,&elasf[i]);
            }
          }
        }
        for(std::vector<MInterfaceElement*>::iterator it=elasf[i].gib.begin(); it!=elasf[i].gib.end();++it){
          MInterfaceElement *ie = *it;
          // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
          MElement *em = ie->getElem(0);
          int edge = ie->getEdgeNumber(0);
          int npts_inter=integrator_on_inter.getIntPoints(ie,&GP);
          for(int j=0;j<npts_inter;j++){
            IPnum key=IPnum(em->getNum(),IPnum::createTypeWithTwoInts(edge,j));
            createAndStoreIP(&key,&elasf[i]);
          }
        }

        for (groupOfElements::elementContainer::const_iterator it = elasf[i].g->begin(); it != elasf[i].g->end(); ++it){
          MElement *e = *it;
          int edge = e->getNumEdges();
          int npts_bulk=integrator_on_bulk.getIntPoints(e,&GP);
          for(int j=0;j<npts_bulk;j++){
            IPnum key=IPnum(e->getNum(),IPnum::createTypeWithTwoInts(edge,j));
            createAndStoreIP(&key,&elasf[i]);
          }
        }
      }
    }
    ~AllIPState(){
      for(std::map<std::pair<long int,int>, IPState*>::iterator it=_mapall.begin();it!=_mapall.end();++it)
        delete (*it).second;
    };
    IPState* getIPstate(IPnum *key){
      return _mapall.find(key->get())->second;
    }
    void nextStep() {state ? state=false : state=true;} // change boolvalue

    void copy(const IPState::whichState ws1, const IPState::whichState ws2){
      for(std::map<std::pair<long int,int>,IPState*>::iterator it=_mapall.begin(); it!=_mapall.end();++it){
        IPState *ips = (*it).second;
        IPVariable *ipv1= ips->getState(ws1);
        IPVariable *ipv2= ips->getState(ws2);
        *ipv2 =*ipv1;
      }
    }
    // Fonction to compute all IP for a given state How to template this function ??
/*    template<class T2> void stressAndDefo(IPState::whichState ws,std::vector<DGelasticField> &elasticFields,
                                            dofManager<double> *pAssembler, T2* LagSpace,
                                            QuadratureBase &Integ_Bulk, QuadratureBase &Integ_Boundary)
 {
  DgC0SolverField<SVector3> SField(pAssembler, LagSpace); //used to interpolate and evaluate gradient
  SVector3 val; // value of a vertex displacement
  IntPt *GP;
  for(int i=0;i<elasticFields.size();i++){
    if(!elasticFields[i].getFormulation()){ // edge gauss point cg/dg (just minus element)
      std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
      std::vector<TensorialTraits<double>::HessType> Hess;
      double uem,vem,uep,vep;
      fullMatrix<double> disp;
      for(std::vector<MInterfaceElement*>::iterator it=elasticFields[i].gi.begin(); it!=elasticFields[i].gi.end();++it){
        MInterfaceElement *ie = *it;
        MElement *e = ie->getElem(0);
        // gauss point
        int npts = Integ_Boundary.getIntPoints(ie,&GP);
        // vector with nodal displacement
        int nbdof = LagSpace->getNumKeys(e);
        int nbFF = e->getNumVertices();
        disp.resize(nbdof,1);
        disp.setAll(0.);
        for(int j=0;j<nbFF;j++){
          SField.getVertexDisplacement(e,val,false,j);
          disp(j,0) = val(0);
          disp(j+nbFF,0) = val(1);
          disp(j+2*nbFF,0) = val(2);
        }
        for(int j=0;j<npts;j++){
          // key of gauss point
          //grad value at gauss point
          ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
          LagSpace->gradfuvw(e,uem,vem,0.,Gradm);
          LagSpace->gradfuvw(ie->getElem(1),uep,vep,0.,Gradp);
          LagSpace->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
          LagSpace->hessfuvw(e,uem,vem,0.,Hess);
          // local basis on element is needed to compute the local basis on interfaceelement (normal)
          LocalBasis lbm,lbp;
          lbm.set(e,Gradm);
          lbp.set(ie->getElem(1),Gradp);
          IPnum key(e->getNum(),IPnum::createTypeWithTwoInts(ie->getEdgeNumber(0),j));
          IPState* ips = this->getIPstate(&key);
          IPVariable* ipv = ips->getState(ws);
          ipv->setLocalBasis(ie,Grads,lbm.gett0(),lbp.gett0());
          ipv->computeStressAndDeformation(elasticFields[i].getMaterialLaw(),nbFF,nbdof,disp,Gradm,Hess);
          // appened method in gradfuvw
          Gradm.clear(); Gradp.clear(); Grads.clear(); Hess.clear();
        }
      }
    }
    else{ //edge gauss point full dg
      std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
      std::vector<TensorialTraits<double>::HessType> Hessm,Hessp;
      double uem,vem,uep,vep;
      fullMatrix<double> dispm,dispp;
      for(std::vector<MInterfaceElement*>::iterator it=elasticFields[i].gi.begin(); it!=elasticFields[i].gi.end();++it){
        MInterfaceElement *ie = *it;
        MElement *em = ie->getElem(0);
        MElement *ep = ie->getElem(1);
        // gauss point
        int npts = Integ_Boundary.getIntPoints(ie,&GP);
        // vector with nodal displacement
        int nbdofm = LagSpace->getNumKeys(em);
        int nbdofp = LagSpace->getNumKeys(ep);
        int nbFFm = em->getNumVertices();
        int nbFFp = ep->getNumVertices();
        dispm.resize(nbdofm,1);
        dispm.setAll(0.);
        dispp.resize(nbdofp,1);
        dispp.setAll(0.);
        for(int j=0;j<nbFFm;j++){
          SField.getVertexDisplacement(em,val,true,j);
          dispm(j,0) = val(0);
          dispm(j+nbFFm,0) = val(1);
          dispm(j+2*nbFFm,0) = val(2);
        }
        for(int j=0;j<nbFFp;j++){
          SField.getVertexDisplacement(ep,val,true,j);
          dispp(j,0) = val(0);
          dispp(j+nbFFp,0) = val(1);
          dispp(j+2*nbFFp,0) = val(2);
        }
        for(int j=0;j<npts;j++){
          // key of gauss point
          //grad value at gauss point
          ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
          LagSpace->gradfuvw(em,uem,vem,0.,Gradm);
          LagSpace->hessfuvw(em,uem,vem,0.,Hessm);
          LagSpace->gradfuvw(ep,uep,vep,0.,Gradp);
          LagSpace->hessfuvw(ep,uep,vep,0.,Hessp);
          LagSpace->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
          // local basis on element is needed to compute the local basis on interfaceelement (normal)
          LocalBasis lbm,lbp;
          lbm.set(em,Gradm);
          lbp.set(ep,Gradp);
          IPnum key(em->getNum(),IPnum::createTypeWithTwoInts(ie->getEdgeNumber(0),j));
          IPState* ips = this->getIPstate(&key);
          IPVariable* ipv = ips->getState(ws);
          ipv->setLocalBasis(ie,Grads,lbm.gett0(),lbp.gett0());
          ipv->computeStressAndDeformation(elasticFields[i].getMaterialLaw(),nbFFm,nbdofm,dispm,Gradm,Hessm);

          // plus elem
          IPnum keyp(ep->getNum(),IPnum::createTypeWithTwoInts(ie->getEdgeNumber(1),j));
          ips = this->getIPstate(&keyp);
          ipv = ips->getState(ws);
          ipv->setLocalBasis(ie,Grads,lbm.gett0(),lbp.gett0());
          ipv->computeStressAndDeformation(elasticFields[i].getMaterialLaw(),nbFFp,nbdofp,dispp,Gradp,Hessp);

          // appened method in gradfuvw
          Grads.clear(); Gradm.clear();Gradp.clear(); Hessm.clear(); Hessp.clear();

        }
      }
    }
    // Virtual interface element
    std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
    std::vector<TensorialTraits<double>::HessType> Hess;
    double uem,vem,uep,vep;
    fullMatrix<double> disp;
    for(std::vector<MInterfaceElement*>::iterator it=elasticFields[i].gib.begin(); it!=elasticFields[i].gib.end();++it){
      MInterfaceElement *ie = *it;
      // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
      MElement *e = ie->getElem(0);
      int edge = ie->getEdgeNumber(0);
      int npts_inter=Integ_Boundary.getIntPoints(ie,&GP);
      // vector with nodal displacement
      int nbdof = LagSpace->getNumKeys(e);
      int nbFF = e->getNumVertices();
      disp.resize(nbdof,1);
      disp.setAll(0.);
      for(int j=0;j<nbFF;j++){
        SField.getVertexDisplacement(e,val,elasticFields[i].getFormulation(),j);
        disp(j,0) = val(0);
        disp(j+nbFF,0) = val(1);
        disp(j+2*nbFF,0) = val(2);
      }
      for(int j=0;j<npts_inter;j++){
        IPnum key=IPnum(e->getNum(),IPnum::createTypeWithTwoInts(edge,j));
        // key of gauss point
        //grad value at gauss point
        ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
        LagSpace->gradfuvw(e,uem,vem,0.,Gradm);
        LagSpace->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
        LagSpace->hessfuvw(e,uem,vem,0.,Hess);
        // local basis on element is needed to compute the local basis on interfaceelement (normal)
        LocalBasis lbm;
        lbm.set(e,Gradm);
        IPState* ips = this->getIPstate(&key);
        IPVariable* ipv = ips->getState(ws);
        ipv->setLocalBasis(ie,Grads,lbm.gett0());
        ipv->computeStressAndDeformation(elasticFields[i].getMaterialLaw(),nbFF,nbdof,disp,Gradm,Hess);
        Gradm.clear(); Grads.clear(); Hess.clear();
      }
    }
    // bulk
    for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it){
          MElement *e = *it;
          int edge = e->getNumEdges();
          int npts_bulk=Integ_Bulk.getIntPoints(e,&GP);
          int nbdof = LagSpace->getNumKeys(e);
          int nbFF = e->getNumVertices();
          disp.resize(nbdof,1);
          disp.setAll(0.);
          for(int j=0;j<nbFF;j++){
            SField.getVertexDisplacement(e,val,elasticFields[i].getFormulation(),j);
            disp(j,0) = val(0);
            disp(j+nbFF,0) = val(1);
            disp(j+2*nbFF,0) = val(2);
          }
          for(int j=0;j<npts_bulk;j++){
            IPnum key=IPnum(e->getNum(),IPnum::createTypeWithTwoInts(edge,j));
            //grad value at gauss point
            LagSpace->gradfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
            LagSpace->hessfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Hess);
        // local basis on element is needed to compute the local basis on interfaceelement (normal)
        IPState* ips = this->getIPstate(&key);
        IPVariable* ipv = ips->getState(ws);
        ipv->setLocalBasis(e,Grads);
        ipv->computeStressAndDeformation(elasticFields[i].getMaterialLaw(),nbFF,nbdof,disp,Grads,Hess);
        // appened method in gradfuvw
        Grads.clear(); Hess.clear();
          }
        }

  }
}*/

};
#endif // IPSTATE_H_
