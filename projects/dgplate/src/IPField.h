//
// C++ Interface: terms
//
// Description: Class to compute Internal point
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
# ifndef _IPFIELD_H_
# define _IPFIELD_H_
#include<vector>
#include"DgC0PlateDof.h"
#include"quadratureRules.h"
#include"IPState.h"
#include "displacementField.h"
#include "elementField.h"
template<class T1, class T2>
class IPField : public elementField {
  protected :
    std::vector<T1>* _efield;
    dofManager<double> *_dm;
    T2* _space;
    QuadratureBase *_intBulk;
    QuadratureBase *_intBound;
    AllIPState *_AIPS;
    displacementField *_ufield; // space field ??

    std::pair<MInterfaceElement*,int> ctp; // use to know the crack tip position;

    // function to compute state depends on element type (template in place of dynamic cast ??)
    void compute1statePlatePlaneStress(IPState::whichState ws, T1* ef);
    double getVMPlatePlaneStress(MElement *ele, const IPState::whichState ws,
                                 const int num, const DGelasticField *elas) const;
    double getSigmaWithOperationPlatePlaneStress(MElement *ele, const IPState::whichState ws, const int num,
                                                     const component::enumcomp cmp, const DGelasticField *elas) const;


    void compute1statePlatePlaneStressWTI(IPState::whichState ws, T1* ef);
    double getVMPlatePlaneStressWTI(MElement *ele, const IPState::whichState ws,
                                 const int num, const DGelasticField *elas, const int pos) const;
    double getStressWithOperationPlatePlaneStressWTI(MElement *ele, const IPState::whichState ws, const int num,
                                                     const component::enumcomp cmp, const DGelasticField *elas, const int pos) const;
    const LocalBasis * getStressReducedWTI(MInterfaceElement *iele,const IPState::whichState ws,
                                  const int num,const DGelasticField *elas, reductionElement &stressTensor,const int pos);

    const void getStressReducedWTI(MInterfaceElement *iele,const IPState::whichState ws,
                                  const int num,const DGelasticField *elas, reductionElement &stressTensor,const int pos, const LocalBasis*[2]);
    void compute1statePlatePlaneStressWF(IPState::whichState ws, T1* ef);
    void setBroken(MInterfaceElement *ie,IPState::whichState ws, const int numminus, const int numplus,
                   const double svm,const SolElementType::eltype elt, const double Gc, const double betaML, const bool tension_){

      std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
      IPVariablePlateOIWF* ipv = dynamic_cast<IPVariablePlateOIWF*>((*vips)[numminus]->getState(ws));
      IPVariablePlateOIWF* ipvp = dynamic_cast<IPVariablePlateOIWF*>((*vips)[numminus]->getState(IPState::previous));
      if(!ipvp->getBroken()){
        Msg::Info("Interface element %d is broken at gauss point %d \n",ie->getNum(), numminus);
        // debugging info
        Msg::Info("minus element = %d plus element = %d\n",ie->getElem(0)->getNum(),ie->getElem(1)->getNum());
        Msg::Info("Position of minus interface\n");
        Msg::Info("x0 = %lf y0 = %lf z0 = %lf x1 = %lf y1 = %lf z1 = %lf",ie->getElem(0)->getVertex(0)->x(),
                  ie->getElem(0)->getVertex(0)->y(),ie->getElem(0)->getVertex(0)->z(),
                  ie->getElem(0)->getVertex(1)->x(),ie->getElem(0)->getVertex(1)->y(),ie->getElem(0)->getVertex(1)->z());
        if(tension_) Msg::Info("tension");
        else Msg::Info("compression");
        ctp.first = ie;
        ctp.second = numminus;
        reductionElement nhatmean, mhatmean;
        const LocalBasis* Lb[3];

        IntPt *GP;
        int npts = _intBound->getIntPoints(ie,&GP);
        // compute n0 and m0
        this->getReductionAndLocalBasis(ie,numminus,npts,elt,ws,nhatmean,mhatmean,Lb);
        std::vector<double> dispm;
        std::vector<double> dispp;
        SVector3 ujump;
        double rjump[3];
        dispm.resize(_space->getNumKeys(ie->getElem(0)));
        dispp.resize(_space->getNumKeys(ie->getElem(1)));
        _ufield->get(ie->getElem(0),dispm);
        _ufield->get(ie->getElem(1),dispp);
        // initial jump
        double uem,uep,vem,vep;
        int nbFFm = ie->getElem(0)->getNumVertices();
        int nbFFp = ie->getElem(1)->getNumVertices();
        ie->getuvOnElem(GP[numminus].pt[0],uem,vem,uep,vep);
        std::vector<TensorialTraits<double>::ValType> Valm;
        _space->fuvw(ie->getElem(0),uem, vem,0., Valm);
        std::vector<TensorialTraits<double>::ValType> Valp;
        _space->fuvw(ie->getElem(1),uep, vep,0., Valp);
        std::vector<TensorialTraits<double>::GradType> Gradm,Gradp;
        _space->gradfuvw(ie->getElem(0),uem,vem,0.,Gradm);
        _space->gradfuvw(ie->getElem(1),uep,vep,0.,Gradp);
        displacementjump(Valm,nbFFm,Valp,nbFFp,dispm,dispp,ujump);
        rotationjump(Lb[2],Gradm,ie->getElem(0)->getNumVertices(),Lb[0],Gradp,ie->getElem(1)->getNumVertices(),Lb[1],dispm,dispp,rjump);
        // Store data
        ipv->setBroken(svm,Gc,betaML,nhatmean,mhatmean,ujump,rjump,Lb[2],tension_);
        ipv = dynamic_cast<IPVariablePlateOIWF*>((*vips)[numplus]->getState(ws));
        ipv->setBroken(svm,Gc,betaML,nhatmean,mhatmean,ujump,rjump,Lb[2],tension_);
      }
    }
    // fracture criteria is based on Camacho & Ortiz Modelling of impact damage (Int J. Solids Structures 1996)
    // here this criteria is evaluated on lower and upper fiber of plate.
    void computeFracture(IPState::whichState ws, T1* ef){
      // For now no fracture at extremities (only MInterfaceElement and no VirtualInterfaceElement)
      IntPt *GP;
      int npts;
      //double svm;
      int msimpm1 = ef->getmsimp()-1;
      reductionElement reducedstressTensor, reducedstressTensorHat, reducedstressTensorMean;
      linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
      double sigmac = mlaw->getSigmac();
      double seff, smax, snor, tau;
      bool ift;
      const LocalBasis* lbb[2];
      for(std::vector<MInterfaceElement*>::iterator it=ef->gi.begin(); it!=ef->gi.end();++it){
        MInterfaceElement *ie = *it;
        npts = _intBound->getIntPoints(ie,&GP);
        // TODO no computation if already broken
        for(int j=0;j<npts;j++){
          // mean value (+ and - element)
          this->getStressReducedWTI(ie,ws,j,ef,reducedstressTensor,0,lbb);
          stressReductionHat(reducedstressTensor,lbb[0],reducedstressTensorMean);
          this->getStressReducedWTI(ie,ws,j+npts,ef,reducedstressTensor,0,lbb);
          stressReductionHat(reducedstressTensor,lbb[0],reducedstressTensorHat);
          for(int k=0;k<2;k++)
            for(int kk=0;kk<2;kk++){
              reducedstressTensorMean(k,kk)+=reducedstressTensorHat(k,kk);
              reducedstressTensorMean(k,kk)*=0.5;
            }
          stressTensor stress(reducedstressTensorMean,lbb[1]);
          snor = stress.getComponent(lbb[1]->getOrthonormalVector(1),lbb[1]->getOrthonormalVector(1));
          tau = stress.getComponent(lbb[1]->getOrthonormalVector(0), lbb[1]->getOrthonormalVector(1));
          // sigma eff (Camacho and Ortiz)
          if(snor>=0.){
            seff = sqrt(snor*snor+mlaw->getBeta()*tau*tau);
            smax = snor;
            ift = true;
          }
          else{
            double temp = fabs(tau)-mlaw->getMu()*fabs(snor);
            if (temp >0)
              seff = sqrt(mlaw->getBeta())*temp;
            else
              seff = 0.;
            smax = tau;
            ift=false;
          }
          if(seff> sigmac){
            this->setBroken(ie,ws,j,j+npts,smax,ef->getSolElemType(),mlaw->getGc(),mlaw->getBeta(),ift);
          } // no computation for last Simpson's point if it already broken
          else{
              reducedstressTensorMean.setAll(0.);
              this->getStressReducedWTI(ie,ws,j,ef,reducedstressTensor,msimpm1,lbb);
              stressReductionHat(reducedstressTensor,lbb[0],reducedstressTensorMean);
              this->getStressReducedWTI(ie,ws,j+npts,ef,reducedstressTensor,msimpm1,lbb);
              stressReductionHat(reducedstressTensor,lbb[0],reducedstressTensorHat);
              for(int k=0;k<2;k++)
                for(int kk=0;kk<2;kk++){
                  reducedstressTensorMean(k,kk)+=reducedstressTensorHat(k,kk);
                  reducedstressTensorMean(k,kk)*=0.5;
              }
              stressTensor stress(reducedstressTensorMean,lbb[1]);
              snor = stress.getComponent(lbb[1]->getOrthonormalVector(1),lbb[1]->getOrthonormalVector(1));
              tau = stress.getComponent(lbb[1]->getOrthonormalVector(0), lbb[1]->getOrthonormalVector(1));
              if(snor>=0.){
                seff = sqrt(snor*snor+mlaw->getBeta()*tau*tau);
                smax = snor;
                ift = true;
              }
              else{
                double temp = fabs(tau)-mlaw->getMu()*fabs(snor);
                if (temp >0)
                  seff = sqrt(mlaw->getBeta())*temp;
                else
                  seff = 0.;
                smax = tau;
                ift = false;
              }
              if(seff>sigmac)
                this->setBroken(ie,ws,j,j+npts,smax,ef->getSolElemType(),mlaw->getGc(),mlaw->getBeta(),ift);

          }
        }
      }
    }
  void getReductionByCohesiveLawAndLocalBasis(const MInterfaceElement *iele, const int gaussnum, const int npts,
                                              reductionElement &nhatmean, reductionElement &mhatmean,
                                              const LocalBasis* lbb[3]){
    // find elasticField
    bool flag=false;
    DGelasticField *ef;
    for(int i=0;i<_efield->size();i++){
      for(std::vector<MInterfaceElement*>::iterator it = (*_efield)[i].gi.begin(); it != (*_efield)[i].gi.end(); ++it){
        MInterfaceElement *ie = *it;
        if(ie==iele){
          flag=true;
          break;
        }
      }
      if(flag) {ef=&(*_efield)[i]; break;}
    }
    // get ipv
    IPVariablePlateOIWF *ipv;
    std::vector<IPState*> *vips;
    IPState *ips;
    vips = _AIPS->getIPstate(iele->getNum());
    ips = (*vips)[gaussnum];
    ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));

    // use the law
    linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
    mlaw->getCohesiveReduction(ipv->getm0(),ipv->getn0(),ipv->getDeltan(),ipv->getDeltanmax(),ipv->getDeltat(),ipv->getDeltatmax(),ipv->getDeltac(),ipv->ifTension(),nhatmean,mhatmean);

    // set localBasis
    lbb[0] = ipv->getLocalBasis();
    lbb[2] = ipv->getLocalBasisOfInterface();
    ips = (*vips)[gaussnum+npts];
    ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));
    lbb[1] = ipv->getLocalBasis();
  }
  void getReductionByCohesiveLaw(const MInterfaceElement *iele, const int gaussnum, const int npts,
                                              reductionElement &nhatmean, reductionElement &mhatmean){
    // find elasticField
    bool flag=false;
    DGelasticField *ef;
    for(int i=0;i<_efield->size();i++){
      for(std::vector<MInterfaceElement*>::iterator it = (*_efield)[i].gi.begin(); it != (*_efield)[i].gi.end(); ++it){
        MInterfaceElement *ie = *it;
        if(ie==iele){
          flag=true;
          break;
        }
      }
      if(flag) {ef=&(*_efield)[i]; break;}
    }
    // get ipv
    IPVariablePlateOIWF *ipv;
    std::vector<IPState*> *vips;
    IPState *ips;
    vips = _AIPS->getIPstate(iele->getNum());
    ips = (*vips)[gaussnum];
    ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));

    // use the law
    linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
    mlaw->getCohesiveReduction(ipv->getm0(),ipv->getn0(),ipv->getDeltan(),ipv->getDeltanmax(),ipv->getDeltat(),ipv->getDeltatmax(),
                               ipv->getDeltac(),ipv->ifTension(),
                               nhatmean,mhatmean);
  }


 public :
  enum ElemValue{mean=10001, max=10002, min=10003}; // enum to select a particular value on element
  IPField(std::vector< T1 >* ef,dofManager<double>* pa,T2* sp,
          QuadratureBase *intb, QuadratureBase *intb1, GModelWithInterface *pmo, displacementField* uf) : _efield(ef), _dm(pa), _space(sp),
                                                                           _intBulk(intb), _intBound(intb1), _ufield(uf),
                                                                           elementField("stress.msh",1000000,1,
                                                                                        elementField::ElementData,true){
  ctp.first = NULL; ctp.second = -1;
  system ("rm crackTipPosition.csv");
  // Creation of storage for IP data
  _AIPS = new AllIPState(pmo, *_efield, *_intBulk, *_intBound);
  // compute the number of element (FIX IT TODO ??)
  long int nelem=0;
  for (unsigned int i = 0; i < _efield->size(); ++i)
    for (groupOfElements::elementContainer::const_iterator it = (*_efield)[i].g->begin(); it != (*_efield)[i].g->end(); ++it)
      nelem++;
  this->setTotElem(nelem);
  this->buildView(*_efield,0.,0,"VonMises",-1,false);
  this->buildView(*_efield,0.,0,"sigmaxx",0,false);
  this->buildView(*_efield,0.,0,"sigmayy",1,false);
  this->buildView(*_efield,0.,0,"tauxy",3,false);
  }
  AllIPState* getAips() const {return _AIPS;}
  ~IPField(){delete _AIPS;}
  void compute1state(IPState::whichState ws){
    for(int i=0;i<_efield->size();i++){
      switch((*_efield)[i].getSolElemType()){
        case SolElementType::PlatePlaneStress :
          this->compute1statePlatePlaneStress(ws,&(*_efield)[i]);
          break;
        case SolElementType::PlatePlaneStressWTI :
        this->compute1statePlatePlaneStressWTI(ws,&(*_efield)[i]);
          break;
        case SolElementType::PlatePlaneStressWF :
          this->compute1statePlatePlaneStressWF(ws,&(*_efield)[i]); // The state is compute in the same way as there is not fracture
      }                                                              // an other function computes fracture
    }

  }
  // On element only ?? Higher level to pass the associated element (pass edge num)
  double getVonMises(MElement *ele, const IPState::whichState ws, const int num, const int pos=0) const{
    double svm;
    // Find elastic field of the element
    bool flag=false;
    DGelasticField *ef;
    for(int i=0;i<_efield->size();i++){
      for(groupOfElements::elementContainer::const_iterator it = (*_efield)[i].g->begin(); it != (*_efield)[i].g->end(); ++it){
        MElement *e = *it;
        if(e==ele){
          flag=true;
          break;
        }
      }
      if(flag) {ef=&(*_efield)[i]; break;}
    }
    // function depends on element type
    switch(ef->getSolElemType()){
      case SolElementType::PlatePlaneStress :
        svm = getVMPlatePlaneStress(ele,ws,num,ef);
        break;
      case SolElementType::PlatePlaneStressWTI :
        svm = getVMPlatePlaneStressWTI(ele,ws,num,ef,0);
        break;
      case SolElementType::PlatePlaneStressWF :
        svm = getVMPlatePlaneStressWTI(ele,ws,num,ef,0);
        break;
      default :
        Msg::Error("Function getVonMises doesn't exist for element type : %d",ef->getSolElemType());
        svm = 0.;
      }
    return svm;
  }

  // get value with a operation
  double getStressWithOperation(MElement *ele, const IPState::whichState ws, const int num, const component::enumcomp cmp, const int pos=0) const{
    double sig;
    // Find elastic field of the element
    bool flag=false;
    DGelasticField *ef;
    for(int i=0;i<_efield->size();i++){
      for(groupOfElements::elementContainer::const_iterator it = (*_efield)[i].g->begin(); it != (*_efield)[i].g->end(); ++it){
        MElement *e = *it;
        if(e==ele){
          flag=true;
          break;
        }
      }
      if(flag) {ef=&(*_efield)[i]; break;}
    }
    // function depends on element type
    switch(ef->getSolElemType()){
      case SolElementType::PlatePlaneStress :
        sig = this->getSigmaWithOperationPlatePlaneStress(ele,ws,num,cmp,ef);
        break;
      case SolElementType::PlatePlaneStressWTI :
        sig = this->getStressWithOperationPlatePlaneStressWTI(ele,ws,num,cmp,ef,0);
        break;
      case SolElementType::PlatePlaneStressWF :
        sig = this->getStressWithOperationPlatePlaneStressWTI(ele,ws,num,cmp,ef,0);
        break;
      default :
        Msg::Error("Function getSigmaWithOperation doesn't exist for element type : %d",ef->getSolElemType());
        sig = 0.;
      }
    return sig;
  }

  // function to archive
  virtual void get(MElement *ele, std::vector<double> &stress, const int cc=-1){
    switch(cc){
      case -1 :
        stress[0]= this->getVonMises(ele,IPState::current,mean,0);
        break;
      case 0 :
        stress[0] = this->getStressWithOperation(ele,IPState::current,mean,component::xx,0);
        break;
      case 1 :
        stress[0] = this->getStressWithOperation(ele,IPState::current,mean,component::yy);
        break;
      case 2 :
        stress[0] = this->getStressWithOperation(ele,IPState::current,mean,component::zz);
        break;
      case 3 :
        stress[0] = this->getStressWithOperation(ele,IPState::current,mean,component::xy);
        break;
    }
  }

  void evalFracture(IPState::whichState ws){
    for(int i=0;i<_efield->size();i++){
      switch((*_efield)[i].getSolElemType()){
        case SolElementType::PlatePlaneStressWF : // other case no fracture
          computeFracture(ws,&(*_efield)[i]);
          break;
      }
    }
  }

  // retrieve a vector with fractured gauss point for an element
  void getBroken(MInterfaceElement* iele, const SolElementType::eltype elt, std::vector<bool> &vectB, std::vector<bool> &vectfB){
    IntPt *GP;
    int npts = _intBound->getIntPoints(iele,&GP);
    vectB.resize(2*npts);
    vectfB.resize(2*npts);
    for(int i=0;i<2*npts;i++) vectfB[i]=false;
    // cG/dG case no fracture
    if(elt != SolElementType::PlatePlaneStressWF)
      for(int i=0;i<2*npts;i++) vectB[i]=false;
    else{
      // retrieve information from AIPS
      IPVariablePlateOIWF *ipv;
      std::vector<IPState*> *vips;
      IPState *ips;
      vips = _AIPS->getIPstate(iele->getNum());
      for(int j=0;j<2*npts;j++){
        ips = (*vips)[j];
        ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));
        vectB[j] = ipv->getBroken();
        if(vectB[j]){
          if(ipv->getDeltan()<0.) vectfB[j] =true;
        }
      }
    }
  }
  bool getBroken(const MInterfaceElement *iele, const int numgauss, const SolElementType::eltype elt, const IPState::whichState ws){
    bool t = false;
    if(elt == SolElementType::PlatePlaneStressWF){
      IPVariablePlateOIWF *ipv;
      std::vector<IPState*> *vips;
      IPState *ips;
      vips = _AIPS->getIPstate(iele->getNum());
      ips = (*vips)[numgauss];
      ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(ws));
      t = ipv->getBroken();
    }
    return t;
  }

  // Interaction with Aips
  void copy(IPState::whichState source, IPState::whichState dest){_AIPS->copy(source,dest);}
  void nextStep(){
    // If fracture deltamax must be updated before update
    // deltamax is set in current state just before it becomes previous state;
    for(int i=0;i<_efield->size();i++){
      if((*_efield)[i].getSolElemType() == SolElementType::PlatePlaneStressWF){
        // loop on interfaceElement
        IPVariablePlateOIWF *ipv;
        std::vector<IPState*> *vips;
        IPState *ips;
        for(std::vector<MInterfaceElement*>::iterator it=(*_efield)[i].gi.begin(); it!=(*_efield)[i].gi.end();++it){
          MInterfaceElement *ie = *it;
          vips = _AIPS->getIPstate(ie->getNum());
          for(int j=0;j<vips->size();j++){
            ips = (*vips)[j];
            ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));
            if(ipv->getBroken()) ipv->updatedeltamax();
          }
        }
      }
    }
    _AIPS->nextStep();
  }

  // initial broken
  void initialBroken(GModelWithInterface* pModel, std::vector<int> &vnumphys);
  void initialBroken(MInterfaceElement *iele, materialLaw *mlaw);

  // reduction element
  void getStressReduction(MElement *ele,const int gaussnum, SolElementType::eltype et,IPState::whichState ws,
                          reductionElement &nalpha);
  void getMomentReduction(MElement *ele,const int gaussnum, SolElementType::eltype et,IPState::whichState ws,
                          reductionElement &malpha);
  void getReduction(MElement *ele,const int gaussnum, SolElementType::eltype et,IPState::whichState ws,
                          reductionElement &nalpha,reductionElement &malpha );

  const LocalBasis* getReductionAndLocalBasis(MElement *ele,const int gaussnum, SolElementType::eltype et,IPState::whichState ws,
                          reductionElement &nalpha,reductionElement &malpha );

  void getStressReduction(MInterfaceElement *iele, const int gaussnum, const int numOfGaussPoint, SolElementType::eltype et,
                                         IPState::whichState ws, reductionElement &nhatmean);
  void getMomentReduction(MInterfaceElement *iele, const int gaussnum, const int numOfGaussPoint, SolElementType::eltype et,
                                         IPState::whichState ws, reductionElement &mhatmean);
  void getMomentReductionAndLocalBasis(MInterfaceElement *iele, const int gaussnum, const int numOfGaussPoint, SolElementType::eltype et,
                                         IPState::whichState ws, reductionElement &mhatmean, const LocalBasis *lb[3]);

  void getReductionAndLocalBasis(MInterfaceElement *iele, const int gaussnum, const int numOfGaussPoint, SolElementType::eltype et,
                                         IPState::whichState ws, reductionElement &nhatmean, reductionElement &mhatmean,
                                         const LocalBasis *lb[3]);
  void getVirtualMomentReductionAndLocalBasis(MInterfaceElement *iele, const int gaussnum, const int numOfGaussPoint,
                                              SolElementType::eltype et, IPState::whichState ws, reductionElement &mhatmean,
                                              const LocalBasis **lb);
  void getVirtualMomentReduction(MInterfaceElement *iele, const int gaussnum, const int numOfGaussPoint,
                                              SolElementType::eltype et, IPState::whichState ws, reductionElement &mhatmean);

  void getStress(const MElement *ele, const int gaussnum, IPState::whichState st, double stress[6]);
  void getStressMembrane(const MElement *ele, const int gaussnum, IPState::whichState st, double stress[6], double & h);
  void getStressBending(const MElement *ele, const int gaussnum, IPState::whichState st, double stress[6], double & h);
  void getStress(const MElement *ele, const int gaussnum, IPState::whichState st, double stressM[6], double stressB[6], double & h);
  const LocalBasis* getStressAndLocalBasis(const MElement *ele, const int gaussnum, IPState::whichState st, double stressM[6], double stressB[6], double & h);
  void getStressAndLocalBasis(const MInterfaceElement *ele, const int gaussnum, IPState::whichState st, double stressM[6], double stressB[6], double & h,
                              const LocalBasis* lb[3]);

  // WTI
  void getStress(const MElement *ele, const int gaussnum, IPState::whichState st, std::vector<tab6> &stress);
  void getStress(const MElement *ele, const int gaussnum, IPState::whichState st, std::vector<tab6> &stress,
                 std::vector<double> &hsimp);
  LocalBasis* getStressAndLocalBasis(const MElement *ele, const int gaussnum, IPState::whichState st, std::vector<tab6> &stress,
                                     std::vector<double> &hsimp);
  void getStressAndLocalBasis(const MInterfaceElement *ele, const int gaussnum, IPState::whichState st, std::vector<tab6> &stress,
                              std::vector<double> &hsimp, const LocalBasis *lb[3]);

  // WF (matrix by numerical perturbation)
  void getReductionFracture(const MInterfaceElement *iele, const int numgauss, const int npts,
                            const SolElementType::eltype elemtype, const std::vector<double> disp, const int nbFF_m, const int nbdofm,
                            const std::vector<TensorialTraits<double>::ValType> &Vals_m,
                            const std::vector<TensorialTraits<double>::GradType> &Grads_m,
                            const int nbFF_p,const std::vector<TensorialTraits<double>::ValType> &Vals_p,
                            const std::vector<TensorialTraits<double>::GradType> &Grads_p,
                            reductionElement &nhatmean, reductionElement &mhatmean) const;

  // Function to get data (must be removed)

  // works only without crack branching
  void archCrackTipPosition(const double time){
    if(ctp.first){ // otherwise no crack
      // Compute the position of crack tip
      MInterfaceElement *ie = ctp.first;
      double x0,y0,z0,x1,y1,z1;
      x0 = ie->getVertex(0)->x();
      x1 = ie->getVertex(1)->x();
      y0 = ie->getVertex(0)->y();
      y1 = ie->getVertex(1)->y();
      z0 = ie->getVertex(0)->z();
      z1 = ie->getVertex(1)->z();
      double length = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
      IntPt *GP;
      int npts = _intBound->getIntPoints(ie,&GP);
      // abscisse on gauss point
      double u = GP[ctp.second].pt[0];
      double ut = (1+u)*length/2.;
      double xct = x0 + abs(x1-x0)/length*ut;
      double yct = y0 + abs(y1-y0)/length*ut;
      double zct = z0 + abs(z1-z0)/length*ut;
      FILE *FP = fopen("crackTipPosition.csv","a");
      fprintf(FP,"%lf;%lf;%lf;%lf\n",time,xct,yct,zct);
      fclose(FP);
    }

  }
/*  void getData(int numelem,int numgauss, materialLaw *mlaw, double &N, double &M, double &du, double &dr,double &deltan){
    // It's for fracture so type elem and material law are known
    IPVariablePlateOIWF *ipv;
    std::vector<IPState*> *vips;
    reductionElement nhatmean;
    reductionElement mhatmean;
    IPState *ips;
    vips = _AIPS->getIPstate(numelem);
    ips = (*vips)[numgauss];
    ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));
    deltan = ipv->getDeltan();
    double deltat = ipv->getDeltat();
    dr = ipv->getDeltarnor() ;
    du = ipv->getDeltaunor();
    linearElasticLawPlaneStressWithFracture *mlaw1 = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(mlaw);
    mlaw1->getCohesiveReduction(ipv->getm0(), ipv->getn0(),deltan, ipv->getDeltanmax(),deltat, ipv->getDeltatmax, ipv->getDeltac(),ipv->ifTension(),nhatmean,mhatmean);
    M = mhatmean(1,1);
    N = nhatmean(1,1);
  }*/
};
#endif // IPField
