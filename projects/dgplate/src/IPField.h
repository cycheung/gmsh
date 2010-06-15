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

    // function to compute state depends on element type (template in place of dynamic cast ??)
    void compute1statePlatePlaneStress(IPState::whichState ws, T1* ef);
    double getVMPlatePlaneStress(MElement *ele, const IPState::whichState ws,
                                 const int num, const DGelasticField *elas) const;

    void compute1statePlatePlaneStressWTI(IPState::whichState ws, T1* ef);
    double getVMPlatePlaneStressWTI(MElement *ele, const IPState::whichState ws,
                                 const int num, const DGelasticField *elas, const int pos) const;
    const LocalBasis * getStressTensorWTI(MInterfaceElement *iele,const IPState::whichState ws,
                                  const int num,const DGelasticField *elas, reductionElement &stressTensor,const int pos);

    void compute1statePlatePlaneStressWF(IPState::whichState ws, T1* ef);
    void setBroken(MInterfaceElement *ie,IPState::whichState ws, const int numminus, const int numplus,
                   const double svm,const SolElementType::eltype elt, const double Gc){

      std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
      IPVariablePlateOIWF* ipv = dynamic_cast<IPVariablePlateOIWF*>((*vips)[numminus]->getState(ws));
      IPVariablePlateOIWF* ipvp = dynamic_cast<IPVariablePlateOIWF*>((*vips)[numminus]->getState(IPState::previous));
      if(!ipvp->getBroken()){
        Msg::Info("Interface element %d is broken\n",ie->getNum());
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
        ipv->setBroken(svm,Gc,nhatmean,mhatmean,ujump,rjump,Lb[2]);
        ipv = dynamic_cast<IPVariablePlateOIWF*>((*vips)[numplus]->getState(ws));
        ipv->setBroken(svm,Gc,nhatmean,mhatmean,ujump,rjump,Lb[2]);
      }
    }
    // fracture criteria is based on VM stress (fracture if s_VM> s_c)
    void computeFracture(IPState::whichState ws, T1* ef){
      // For now no fracture at extremities (only MInterfaceElement and no VirtualInterfaceElement)
      IntPt *GP;
      int npts;
      //double svm;
      double sc = ef->getSigmaC();
      int msimpm1 = ef->getmsimp()-1;
      reductionElement stressTensor, stressTensorHat,stressTensorMean;

      const LocalBasis *lb;
      for(std::vector<MInterfaceElement*>::iterator it=ef->gi.begin(); it!=ef->gi.end();++it){
        MInterfaceElement *ie = *it;
        npts = _intBound->getIntPoints(ie,&GP);
        // TODO no computation if already broken
        for(int j=0;j<npts;j++){
          // mean value (+ and - element)
          //svm = this->getVMPlatePlaneStressWTI(ie,ws,j,ef,0);
          //svm+= this->getVMPlatePlaneStressWTI(ie,ws,j+npts,ef,0);
          //svm*=0.5;
          lb = this->getStressTensorWTI(ie,ws,j,ef,stressTensor,0);
          stressReductionHat(stressTensor,lb,stressTensorMean);
          lb = this->getStressTensorWTI(ie,ws,j+npts,ef,stressTensor,0);
          stressReductionHat(stressTensor,lb,stressTensorHat);
          for(int k=0;k<2;k++)
            for(int kk=0;kk<2;kk++){
              stressTensorMean(k,kk)+=stressTensorHat(k,kk);
              stressTensorMean(k,kk)*=0.5;
            }

          // if(svm>sc) // doesn't work fracture in compression ??
          if(stressTensorMean(1,1)>sc){ // Fracture in mode I only
            linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
            this->setBroken(ie,ws,j,j+npts,stressTensorMean(1,1),ef->getSolElemType(),mlaw->getGc());
          } // no computation for last Simpson's point if it already broken
          else{
//            svm = this->getVMPlatePlaneStressWTI(ie,ws,j,ef,msimpm1);
//            svm+= this->getVMPlatePlaneStressWTI(ie,ws,j+npts,ef,msimpm1);
//            svm*=0.5;
              stressTensorMean.setAll(0.);
              lb = this->getStressTensorWTI(ie,ws,j,ef,stressTensor,msimpm1);
              stressReductionHat(stressTensor,lb,stressTensorMean);
              lb = this->getStressTensorWTI(ie,ws,j+npts,ef,stressTensor,msimpm1);
              stressReductionHat(stressTensor,lb,stressTensorHat);
              for(int k=0;k<2;k++)
                for(int kk=0;kk<2;kk++){
                  stressTensorMean(k,kk)+=stressTensorHat(k,kk);
                  stressTensorMean(k,kk)*=0.5;
                }
            //if(svm>sc){ // doesn't work fracture in compression
            if(stressTensorMean(1,1)>sc){ // fracture in mode I
              //Msg::Info("Interface element %d is upper broken\n",ie->getNum());
              linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
              this->setBroken(ie,ws,j,j+npts,stressTensorMean(1,1),ef->getSolElemType(),mlaw->getGc());
            }
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
    mlaw->getCohesiveReduction(ipv->getM0(),ipv->getN0(),ipv->getDelta(),ipv->getDeltamax(),ipv->getDeltac(),nhatmean,mhatmean);

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
    mlaw->getCohesiveReduction(ipv->getM0(),ipv->getN0(),ipv->getDelta(),ipv->getDeltamax(),ipv->getDeltac(),nhatmean,mhatmean);
  }


 public :
  enum ElemValue{mean=10001, max=10002, min=10003}; // enum to select a particular value on element
  IPField(std::vector< T1 >* ef,dofManager<double>* pa,T2* sp,
          QuadratureBase *intb, QuadratureBase *intb1, GModelWithInterface *pmo, displacementField* uf) : _efield(ef), _dm(pa), _space(sp),
                                                                           _intBulk(intb), _intBound(intb1), _ufield(uf),
                                                                           elementField("stressVM.msh",1000000,1,
                                                                                        elementField::ElementData,"VonMises",true){
  // Creation of storage for IP data
  _AIPS = new AllIPState(pmo, *_efield, *_intBulk, *_intBound);
  // compute the number of element (FIX IT TODO ??)
  long int nelem=0;
  for (unsigned int i = 0; i < _efield->size(); ++i)
    for (groupOfElements::elementContainer::const_iterator it = (*_efield)[i].g->begin(); it != (*_efield)[i].g->end(); ++it)
      nelem++;
  this->setTotElem(nelem);
  this->buildView(*_efield,0.,0,false);
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

  // function to archive
  virtual void get(MElement *ele, std::vector<double> &stress, const int cc=-1){
    stress[0]= this->getVonMises(ele,IPState::current,max,0);
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
  void getBroken(MInterfaceElement* iele, SolElementType::eltype elt, std::vector<bool> &vectB, std::vector<bool> &vectfB){
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
          double delta = ipv->getDelta();
          double deltac = ipv->getDeltac();
          if(delta >= deltac) vectfB[j] =true;
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
                                              const LocalBasis *lb[3]);

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
  void getData(int numelem,int numgauss, materialLaw *mlaw, double &N, double &M, double &du, double &dr,double &delta){
    // It's for fracture so type elem and material law are known
    IPVariablePlateOIWF *ipv;
    std::vector<IPState*> *vips;
    reductionElement nhatmean;
    reductionElement mhatmean;
    IPState *ips;
    vips = _AIPS->getIPstate(numelem);
    ips = (*vips)[numgauss];
    ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));
    delta = ipv->getDelta();
    dr = ipv->getDeltar();
    du = ipv->getDeltaunor();
    linearElasticLawPlaneStressWithFracture *mlaw1 = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(mlaw);
    mlaw1->getCohesiveReduction(ipv->getM0(), ipv->getN0(),delta, ipv->getDeltamax(),ipv->getDeltac(),nhatmean,mhatmean);
    M = mhatmean(1,1);
    N = nhatmean(1,1);
  }
};
#endif // IPField
