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
    // How to remove this ? (Polymorphic (dynamic cast impossible if no virtual functions)).
    virtual void bidon(){};

};
class IPVariablePlate : public IPVariable{
  protected :
    tab6 sigmaMembrane;
    tab6 sigmaBending;
    tab6 epsilon;
    tab6 rho;
    LocalBasis lb;
    double _h;
  public :
    double getThickness() const {return _h;}
    void computeStressAndDeformation(linearElasticLawPlaneStress*,const int, const int, const std::vector<double>&disp,
                                     const std::vector<TensorialTraits<double>::GradType>&,
                                     const std::vector<TensorialTraits<double>::HessType>&);
    // On interface the local basis of Element is used to compute and it's given as argument
    void computeStressAndDeformation(linearElasticLawPlaneStress*,const LocalBasis*,const int, const int,
                                     const std::vector<double>&disp, const std::vector<TensorialTraits<double>::GradType>&,
                                     const std::vector<TensorialTraits<double>::HessType>&);
    // The constructor doesn't take any argument (A set function allow to compute initial data)
    IPVariablePlate(double h) : IPVariable(),_h(h) {
      for(int i=0;i<6;i++){
        sigmaMembrane[i]=0.;
        sigmaBending[i]=0.;
        epsilon[i]=0.;
        rho[i]=0.;
      }
    }
    ~IPVariablePlate(){}
    // copie constructor
    IPVariablePlate(const IPVariablePlate &source):IPVariable(source)
    {
      for(int i=0;i<6;i++){
        sigmaMembrane[i]=source.getSigmaMembrane(i);
        sigmaBending[i]=source.getSigmaBending(i);
        epsilon[i]=source.getEpsMembrane(i);
        rho[i]=source.getRho(i);
      }
      lb=source.lb;
      _h=source._h;
    }
    IPVariablePlate &operator = (const IPVariablePlate &source){
      IPVariable::operator = (source);
      for(int i=0;i<6;i++){
        sigmaMembrane[i]=source.getSigmaMembrane(i);
        sigmaBending[i]=source.getSigmaBending(i);
        epsilon[i]=source.getEpsMembrane(i);
        rho[i]=source.getRho(i);
      }
      lb=source.lb;
      _h=source._h;
      return *this;
    }
    virtual double getSigma(const int i) const{return sigmaMembrane[i]+sigmaBending[i];}
    virtual double getSigmaMembrane(const int i) const{return sigmaMembrane[i];}
    virtual void getSigmaMembrane(double sig[6]) const{for(int i=0;i<6;i++) sig[i]=sigmaMembrane[i]; }
    virtual void getSigmaBending(double sig[6]) const{for(int i=0;i<6;i++) sig[i]=sigmaBending[i]; }
    virtual double getSigmaBending(const int i) const{return sigmaBending[i];}
    virtual double getEpsilon(const int i) const{return epsilon[i]+rho[i];}
    virtual double getEpsMembrane(const int i) const{return epsilon[i];}
    virtual double getRho(const int i) const{return rho[i];}
    virtual void getSigma(double sig[6]){for(int i=0;i<6;i++) sig[i]=sigmaMembrane[i]+sigmaBending[i];}
    virtual void getEpsilon(double eps[6]){for(int i=0;i<6;i++) eps[i]=epsilon[i]+rho[i];}
    virtual void setLocalBasis(MElement *e, const std::vector<TensorialTraits<double>::GradType> &Grads){lb.set(e,Grads);}
    virtual void setLocalBasis(const LocalBasis &lbe){lb = lbe;}
    virtual const LocalBasis * getLocalBasis() const {return &lb;}
};

// Idem on interface Element (the localbasis of interface must be store)
class IPVariablePlateOnInterface : public IPVariablePlate{
 protected :
   LocalBasis lbs;
 public :
   IPVariablePlateOnInterface(const double th) : IPVariablePlate(th){}
   ~IPVariablePlateOnInterface(){};
    IPVariablePlateOnInterface(const IPVariablePlateOnInterface &source):IPVariablePlate(source)
    {
      lbs = source.lbs;
    }
    IPVariablePlateOnInterface &operator = (const IPVariablePlateOnInterface &source){
      IPVariablePlate::operator = (source);
      lbs = source.lbs;
      return *this;
    }
    virtual const LocalBasis* getLocalBasisOfInterface() const{return &lbs;}
    virtual const LocalBasis* getLocalBasis() const{return &lb;}
    void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads,
                               const SVector3 &t0p, const SVector3 &t0m)
    {
      lbs.set(ie,Grads,t0p,t0m);
    }
    virtual void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads,
                               const SVector3 &t0m){
      lbs.set(ie,Grads,t0m);
    }
    virtual void setLocalBasisOfInterface(const LocalBasis &lba){lbs = lba;}
    virtual void setLocalBasis(const LocalBasis &lbe){lb = lbe;}
};

class IPVariablePlateWithThicknessIntegration : public IPVariable{
  protected :
    std::vector<tab6> sigma;
    std::vector<tab6> epsilon;
    //std::vector<LocalBasis> lb;
    // Same LocalBasis for all Simpson's point ??
    LocalBasis lb;
    std::vector<double> zsimp;
    short int nsimp; // Give the number of Simpson's points avoid vect.size()
  public :
    // get operation
    short int getNumSimp() const{return nsimp;}
    double getThickness() const {return zsimp[nsimp-1]-zsimp[0];}
    double getSigma(const int i,const int j) const {return sigma[i][j];}
    double getEpsilon(const int i,const int j)const {return epsilon[i][j];}
    double getSigmaNeutral(const int i) const {return sigma[nsimp/2][i];}
    double getEpsilonNeutral(const int i) const{return epsilon[nsimp/2][i];}
    double getSigmaUpper(const int i) const{return sigma[nsimp-1][i];}
    double getEpsilonUpper(const int i) const{return epsilon[nsimp-1][i];}
    double getSigmaLower(const int i) const{return sigma[0][i];}
    double getEpsilonLower(const int i) const{return epsilon[0][i];}
    virtual void getSigma(std::vector<tab6> &sig){
      for(int j=0;j<nsimp;j++)
        for(int i=0;i<6;i++) sig[j][i]=sigma[j][i];
    }
    virtual void getEpsilon(std::vector<tab6> &eps){
      for(int j=0;j<nsimp;j++)
        for(int i=0;i<6;i++) eps[j][i]=epsilon[j][i];
    }
    virtual void getSimpsonPoint(std::vector<double> & zsi) const{for(int i=0;i<nsimp;i++) zsi[i]=zsimp[i];}
    // Constructor
    IPVariablePlateWithThicknessIntegration(const short int numsimp, const double h) : nsimp(numsimp){
      sigma.resize(numsimp);
      epsilon.resize(numsimp);
      zsimp.resize(numsimp);
      // Compute of position whichi is f(numpsimp,h)
      double hsimp = double(h)/double(numsimp-1); // numsimp is always odd by build of elasticfield
      double halfh = h/2.;
      for(int i=0;i<numsimp;i++)
        zsimp[i] = -halfh + i*hsimp;
    }

    ~IPVariablePlateWithThicknessIntegration(){}
    IPVariablePlateWithThicknessIntegration& operator=(const IPVariablePlateWithThicknessIntegration &source)
    {
      nsimp = source.getNumSimp();
      sigma.resize(nsimp);
      epsilon.resize(nsimp);
      zsimp.resize(nsimp);
      lb = source.lb;
      for(int i=0; i<nsimp;i++){
        sigma[i] = source.sigma[i];
        epsilon[i] = source.epsilon[i];
        zsimp[i] = source.zsimp[i];
      }
    }
    IPVariablePlateWithThicknessIntegration(const IPVariablePlateWithThicknessIntegration& source){
      nsimp = source.getNumSimp();
      sigma.resize(nsimp);
      epsilon.resize(nsimp);
      zsimp.resize(nsimp);
      lb = source.lb;
      for(int i=0; i<nsimp;i++){
        sigma[i] = source.sigma[i];
        epsilon[i] = source.epsilon[i];
        zsimp[i] = source.zsimp[i];
      }
    }
    // operation on IP variables.
    void computeStressAndDeformation(linearElasticLawPlaneStress*,const int, const int, const std::vector<double>&disp,
                                     const std::vector<TensorialTraits<double>::GradType>&,
                                     const std::vector<TensorialTraits<double>::HessType>&);
    // On interface the local basis of Element is used to compute and it's given as argument
    void computeStressAndDeformation(linearElasticLawPlaneStress*,const LocalBasis*,const int, const int,
                                     const std::vector<double>&disp, const std::vector<TensorialTraits<double>::GradType>&,
                                     const std::vector<TensorialTraits<double>::HessType>&);
    virtual void setLocalBasis(MElement *e, const std::vector<TensorialTraits<double>::GradType> &Grads){
      lb.set(e,Grads);
    }
/*    virtual void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads, const SVector3 &t0p, const SVector3 &t0m)
    {
      lb.set(ie,Grads,t0p,t0m);
    }
    virtual void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads, const SVector3 &t0m){
      lb.set(ie,Grads,t0m);
    }
    virtual void setlocalBasis(const LocalBasis *lbe){
      lb = *lbe;
    }*/
    LocalBasis * getLocalBasis(){return &lb;}
};

class IPVariablePlateWithThicknessIntegrationOI : public IPVariablePlateWithThicknessIntegration{
 protected :
   LocalBasis lbs;
 public :
   IPVariablePlateWithThicknessIntegrationOI(const short int numsimp,
                                             const double h) : IPVariablePlateWithThicknessIntegration(numsimp,h){}
   ~IPVariablePlateWithThicknessIntegrationOI(){}
    IPVariablePlateWithThicknessIntegrationOI(const IPVariablePlateWithThicknessIntegrationOI &source):
                                              IPVariablePlateWithThicknessIntegration(source)
    {
      lbs = source.lbs;
    }
    IPVariablePlateWithThicknessIntegrationOI &operator = (const IPVariablePlateWithThicknessIntegrationOI &source){
      IPVariablePlateWithThicknessIntegration::operator = (source);
      lbs = source.lbs;
      return *this;
    }
    void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads,
                               const SVector3 &t0p, const SVector3 &t0m)
    {
      lbs.set(ie,Grads,t0p,t0m);
    }
    void setLocalBasis(MInterfaceElement *ie,const std::vector<TensorialTraits<double>::GradType> &Grads,
                               const SVector3 &t0m){
      lbs.set(ie,Grads,t0m);
    }
    void setLocalBasis(const LocalBasis &lba){lb = lba;}
    void setLocalBasisOfInterface(const LocalBasis &lba){lbs = lba;}
    virtual const LocalBasis* getLocalBasisOfInterface() const{return &lbs;}
};

class IPVariablePlateOIWF : public IPVariablePlateWithThicknessIntegrationOI{
 protected :
  // Add some variables use in case of fracture
  std::vector<SVector3> n0;
  std::vector<SVector3> m0;
  bool broken;
  double ujump, rjump, ujump0, rjump0, delta, deltac, delta0, deltamax, sigmac, beta, _M0, _N0; // _N0 usefull ??

 public :
   IPVariablePlateOIWF(const short int numsimp,
                                             const double h) : IPVariablePlateWithThicknessIntegrationOI(numsimp,h),
                                             broken(false), delta(0.), deltamax(0.){} // no initial fracture (for initial fracture use ipfield)
   ~IPVariablePlateOIWF(){}
    IPVariablePlateOIWF(const IPVariablePlateOIWF &source):
                                              IPVariablePlateWithThicknessIntegrationOI(source)
    {
      for(int i=0;i<n0.size();i++)
        n0[i]=source.n0[i];
      for(int i=0;i<m0.size();i++)
        m0[i]=source.m0[i];
      broken = source.broken;
      ujump = source.ujump;
      ujump0= source.ujump0;
      rjump = source.rjump;
      rjump0=source.rjump0;
      delta = source.delta;
      deltac = source.deltac;
      delta0 = source.delta0;
      deltamax = source.deltamax;
      sigmac = source.sigmac;
      beta = source.beta;
    }
    IPVariablePlateOIWF &operator = (const IPVariablePlateOIWF &source){
      IPVariablePlateWithThicknessIntegrationOI::operator = (source);
      for(int i=0;i<n0.size();i++)
        n0[i]=source.n0[i];
      for(int i=0;i<m0.size();i++)
        m0[i]=source.m0[i];
      broken = source.broken;
      ujump = source.ujump;
      ujump0= source.ujump0;
      rjump = source.rjump;
      rjump0=source.rjump0;
      delta = source.delta;
      deltac = source.deltac;
      delta0 = source.delta0;
      deltamax = source.deltamax;
      sigmac = source.sigmac;
      beta = source.beta;
      return *this;
    }
    void setBroken(const double svm, const double Gc, const std::vector<SVector3> nalpha, const std::vector<SVector3> malpha,
                   const double du, const double dr[3])
    {
      if(!broken){ // initialize broken at this gauss point
        double hdiv6 = this->getThickness()/6.;
        n0.resize(nalpha.size());
        m0.resize(malpha.size());
        for(int i=0;i<nalpha.size();i++){
          n0[i][0] = nalpha[i][0] ; n0[i][1] = nalpha[i][1] ; n0[i][2] = nalpha[i][2];
        }
        for(int i=0;i<malpha.size();i++){
          m0[i][0] = malpha[i][0] ; m0[i][1] = malpha[i][1] ; m0[i][2] = malpha[i][2];
        }
        _M0 = malpha[1][1];
        _N0 = nalpha[1][1];
        ujump0 = du;
        ujump =0.;
        rjump0 = dr[1];
        rjump = 0.;
        sigmac = svm;
        delta = 0.;
        deltac = 2*Gc/sigmac;
        deltamax =0.;
        double m0abs;
        if(_M0>0) {m0abs = _M0; rjump0 = -rjump0;}
        else m0abs = -_M0;
        beta = hdiv6*m0abs/(hdiv6*m0abs+_N0);
        delta0 = beta*hdiv6*rjump0 + (1-beta)*ujump0;
        broken = true;
      }
    }
    void updatedeltamax(){if(delta>deltamax) deltamax = delta;}

    void setFracture(const IPVariablePlateOIWF *ipvprev,const double ujump_, const double rjump_[3]){
      int sizen = ipvprev->n0.size();
      int sizem = ipvprev->m0.size();
      n0.resize(sizen);
      m0.resize(sizem);
      for(int i=0; i<sizen;i++) for(int j=0;j<3;j++) {n0[i][j]=ipvprev->n0[i][j];}
      for(int i=0; i<sizem;i++) for(int j=0;j<3;j++) {m0[i][j]=ipvprev->m0[i][j];}
      ujump0 = ipvprev->ujump0;
      rjump0 = ipvprev->rjump0;
      sigmac = ipvprev->sigmac;
      deltac = ipvprev->deltac;
      delta0 = ipvprev->delta0;
      deltamax = ipvprev->deltamax;
      beta = ipvprev->beta;
      broken = ipvprev->broken;
      _M0 = ipvprev->_M0;
      _N0 = ipvprev->_N0;
      ujump = ujump_ - ujump0;
      // tempory set rjump with respect of max(malpha)
      if(_M0>0) rjump = -rjump_[1] - rjump0;
      else rjump = rjump_[1] - rjump0;
      double hdiv6 = this->getThickness()/6.;
      delta = beta*hdiv6*rjump + (1-beta)*ujump;
    }
    double computeDelta(const SVector3 ujump_, const double rjump_[3]) const {
      if(_M0>0) return (1-beta)*(ujump_[0]-ujump0) + beta*this->getThickness()/6.*(-rjump_[1]-rjump0);
      else return (1-beta)*(ujump_[0]-ujump0) + beta*this->getThickness()/6.*(rjump_[1]-rjump0);
    }
    bool getBroken() const{return broken;}
    double getDelta() const{return delta;}
    double getDeltac() const{return deltac;}
    double getM0() const{return _M0;}
    double getN0() const{return _N0;}
    double getDeltamax() const{return deltamax;}
    double getDeltar() const{return rjump;}
};

class IPState{
  protected :
    IPVariable *_initial;     // initial state t=0
    IPVariable *_previous;    // previous step t=n-1
    IPVariable *_current;     // current step
    // pointer on a bool value to choice what vector is current and what vector is previous
    bool *_st;
    SolElementType::eltype _we;
    // pointer to the elasticField of IP
  public :
    enum whichState{initial, previous, current}; //protected enum ??
    IPState(bool *st,const DGelasticField* elasf, bool inter) : _st(st), _we(elasf->getSolElemType())
    {
      switch(_we){
        case SolElementType::PlatePlaneStress :
          if(!inter){
            _initial = new IPVariablePlate(elasf->_h);
            _previous = new IPVariablePlate(elasf->_h);
            _current = new IPVariablePlate(elasf->_h);
          }
          else{
            _initial = new IPVariablePlateOnInterface(elasf->_h);
            _previous = new IPVariablePlateOnInterface(elasf->_h);
            _current = new IPVariablePlateOnInterface(elasf->_h);
          }
          break;
        case SolElementType::PlatePlaneStressWTI :
          if(!inter){
            _initial = new IPVariablePlateWithThicknessIntegration(elasf->getmsimp(), elasf->_h);
            _previous = new IPVariablePlateWithThicknessIntegration(elasf->getmsimp(), elasf->_h);
            _current = new IPVariablePlateWithThicknessIntegration(elasf->getmsimp(), elasf->_h);
          }
          else{
            _initial = new IPVariablePlateWithThicknessIntegrationOI(elasf->getmsimp(), elasf->_h);
            _previous = new IPVariablePlateWithThicknessIntegrationOI(elasf->getmsimp(), elasf->_h);
            _current = new IPVariablePlateWithThicknessIntegrationOI(elasf->getmsimp(), elasf->_h);
          }
          break;
        case SolElementType::PlatePlaneStressWF : // TODO ?? Allow to not integrate on bulk element
          if(!inter){
            _initial = new IPVariablePlateWithThicknessIntegration(elasf->getmsimp(), elasf->_h);
            _previous = new IPVariablePlateWithThicknessIntegration(elasf->getmsimp(), elasf->_h);
            _current = new IPVariablePlateWithThicknessIntegration(elasf->getmsimp(), elasf->_h);
          }
          else{
            _initial = new IPVariablePlateOIWF(elasf->getmsimp(), elasf->_h);
            _previous = new IPVariablePlateOIWF(elasf->getmsimp(), elasf->_h);
            _current = new IPVariablePlateOIWF(elasf->getmsimp(), elasf->_h);
          }
          break;
        default : Msg::Error("Impossible to create IPVariable for this element type\n");
      }
    }
    ~IPState(){
        if(_initial) delete _initial;
        if(_current) delete _current;
        if(_previous) delete _previous;

    }
    IPState(IPState &source) : _we(source._we)
    {
      if(source._we == SolElementType::PlatePlaneStress ){
        IPVariablePlate *psinit = dynamic_cast<IPVariablePlate*>(source.getState(IPState::initial));
        IPVariablePlate *pspre = dynamic_cast<IPVariablePlate*>(source.getState(IPState::previous));
        IPVariablePlate *pscur = dynamic_cast<IPVariablePlate*>(source.getState(IPState::current));
        _initial = new IPVariablePlate(psinit->getThickness());
        _previous = new IPVariablePlate(pspre->getThickness());
        _current = new IPVariablePlate(pscur->getThickness());
      }
      else if((source._we == SolElementType::PlatePlaneStressWTI) or (source._we == SolElementType::PlatePlaneStressWF) ){
        IPVariablePlateWithThicknessIntegration *psinit = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(source.getState(IPState::initial));
        IPVariablePlateWithThicknessIntegration *pspre = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(source.getState(IPState::previous));
        IPVariablePlateWithThicknessIntegration *pscur = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(source.getState(IPState::current));
        _initial = new IPVariablePlateWithThicknessIntegration(psinit->getNumSimp(),psinit->getThickness());
        _previous = new IPVariablePlateWithThicknessIntegration(pspre->getNumSimp(), pspre->getThickness());
        _current = new IPVariablePlateWithThicknessIntegration(pscur->getNumSimp(),pscur->getThickness() );
      }
      _initial->operator = (*(source._initial));
      _previous = source._previous;
      _current = source._current;
      _st = source._st;
    }
    IPState & operator = (const IPState &source){
      _initial->operator = (*(source._initial));
      _initial = source._initial;
      _previous = source._previous;
      _current = source._current;
      _st = source._st;
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
};

// Class to access to the IPState of all gauss point
class AllIPState{
  protected :
    //std::vector < IPState > allIS;
    typedef std::map<long int, std::vector<IPState*> > ipstateContainer;
    ipstateContainer _mapall;
    // flag to switch previous and current (change the pointer in place of copy all variables)
    bool state;
    void createAndStoreIP(const MElement *ele,const int npts, const DGelasticField* pelasf){
      // vector with pointer on IPState object linked to the element
      std::vector<IPState*> tp(npts);
      for(int i=0;i<npts;i++){
        IPState*ips = new IPState(&state,pelasf,false);
	    tp[i]=ips;
      }
      _mapall.insert(std::pair<long int,std::vector<IPState*> >(ele->getNum(),tp));
    }
    void createAndStoreIP(const MInterfaceElement *iele,const int npts, const DGelasticField* pelasf){
      // vector with pointer on IPState object linked to the element
      std::vector<IPState*> tp(npts);
      for(int i=0;i<npts;i++){
        IPState*ips = new IPState(&state,pelasf,true);
	    tp[i]=ips;
      }
      _mapall.insert(std::pair<long int,std::vector<IPState*> >(iele->getNum(),tp));
    }

  public :
    AllIPState(GModelWithInterface *pModel, std::vector<DGelasticField> &elasf, QuadratureBase &integrator_on_bulk, QuadratureBase &integrator_on_inter){
      state = true; // at creation of object state is true
      IntPt *GP; // needed to know the number of gauss point
      for(int i=0;i<elasf.size();i++){
          // loop
          for(std::vector<MInterfaceElement*>::iterator it=elasf[i].gi.begin(); it!=elasf[i].gi.end();++it){
            MInterfaceElement *ie = *it;
     	    // 2* because IP is duplicated for fullDg formulation
            int npts_inter=2*integrator_on_inter.getIntPoints(ie,&GP);
            createAndStoreIP(ie,npts_inter,&elasf[i]);
          }
        // Virtual interface element (no duplication)
        for(std::vector<MInterfaceElement*>::iterator it=elasf[i].gib.begin(); it!=elasf[i].gib.end();++it){
          MInterfaceElement *ie = *it;
          int npts_inter=integrator_on_inter.getIntPoints(ie,&GP);
          createAndStoreIP(ie,npts_inter,&elasf[i]);
        }
        // bulk element
        for (groupOfElements::elementContainer::const_iterator it = elasf[i].g->begin(); it != elasf[i].g->end(); ++it){
          MElement *e = *it;
          int npts_bulk=integrator_on_bulk.getIntPoints(e,&GP);
          createAndStoreIP(e,npts_bulk,&elasf[i]);
        }
      }
    }
    ~AllIPState(){
      for(std::map<long int, std::vector<IPState*> >::iterator it=_mapall.begin();it!=_mapall.end();++it){
        std::vector<IPState*> vips = (*it).second;
	    for(int i=0;i<vips.size();i++)
	      delete vips[i];
      }
    };
    std::vector<IPState*>* getIPstate(const long int num){
      return &(_mapall.find(num)->second);
    }
    void nextStep() {state ? state=false : state=true;} // change boolvalue

    void copy(const IPState::whichState ws1, const IPState::whichState ws2){
      for(std::map<long int,std::vector<IPState*> >::iterator it=_mapall.begin(); it!=_mapall.end();++it){
        std::vector<IPState*> *vips = &((*it).second);
	    for(int i=0;i<vips->size();i++){
          IPVariable *ipv1= (*vips)[i]->getState(ws1);
          IPVariable *ipv2= (*vips)[i]->getState(ws2);
          *ipv2 =*ipv1;
	    }
      }
    }
};
#endif // IPSTATE_H_
