//
// C++ Interface: terms
//
// Description: Definition of function defined in class IPField for ShellPlaneStress element
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "IPField.h"
// Functions to update 1 State
void IPField::compute1stateShellPlaneStress(IPState::whichState ws, partDomain* ef){
  SVector3 val; // value of a vertex displacement
  IntPt *GP;
  std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
  std::vector<TensorialTraits<double>::HessType> Hessm,Hessp;
  double uem,vem,uep,vep;
  std::vector<double> dispm,dispp;
  dgPartDomain* dgdom = dynamic_cast<dgPartDomain*>(ef);
  for(std::vector<MInterfaceElement*>::iterator it=dgdom->gi.begin(); it!=dgdom->gi.end();++it){
    MInterfaceElement *ie = *it;
    MElement *em = ie->getElem(0);
    MElement *ep = ie->getElem(1);
    // gauss point
    int npts = dgdom->getInterfaceGaussIntegrationRule()->getIntPoints(ie,&GP);
    // vector with nodal displacement
    int nbdofm = _space->getNumKeys(em);
    int nbdofp = _space->getNumKeys(ep);
    int nbFFm = em->getNumVertices();
    int nbFFp = ep->getNumVertices();
    _ufield->get(em,dispm);
    _ufield->get(ep,dispp);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
    for(int j=0;j<npts;j++){
      // key of gauss point grad value at gauss point
      ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
      _space->gradfuvw(em,uem,vem,0.,Gradm);
      _space->hessfuvw(em,uem,vem,0.,Hessm);
      _space->gradfuvw(ep,uep,vep,0.,Gradp);
      _space->hessfuvw(ep,uep,vep,0.,Hessp);
      _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      LocalBasis lbm,lbp,lbs;
      lbm.set(em,Gradm,Hessm);
      lbp.set(ep,Gradp,Hessp);
      lbs.set(ie,Grads,lbm.gett0(),lbp.gett0());
      lbm.set_pushForward(&lbs);
      lbp.set_pushForward(&lbs);
      IPState* ips = (*vips)[j];
      IPVariablePlateOnInterface *ipv = dynamic_cast<IPVariablePlateOnInterface*>(ips->getState(ws));
      ipv->setLocalBasis(lbm);
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,nbFFm,nbdofm,dispm,Gradm,Hessm);

      // plus elem
      ips = (*vips)[npts+j];
      ipv = dynamic_cast<IPVariablePlateOnInterface*>(ips->getState(ws));
      ipv->setLocalBasis(lbp) ;
      ipv->setLocalBasisOfInterface(lbs);
      ipv->computeStressAndDeformation(mlaw,nbFFp,nbdofp,dispp,Gradp,Hessp);

      // appened method in gradfuvw
       Grads.clear(); Gradm.clear();Gradp.clear(); Hessm.clear(); Hessp.clear();
    }
    dispm.clear(); dispp.clear();
  }

  // Virtual interface element
  std::vector<double> disp;
  for(std::vector<MInterfaceElement*>::iterator it=dgdom->gib.begin(); it!=dgdom->gib.end();++it){
    MInterfaceElement *ie = *it;
    // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
    MElement *e = ie->getElem(0);
    int edge = ie->getEdgeNumber(0);
    int npts_inter= dgdom->getInterfaceGaussIntegrationRule()->getIntPoints(ie,&GP);
    // vector with nodal displacement
    int nbdof = _space->getNumKeys(e);
    int nbFF = e->getNumVertices();
    _ufield->get(e,disp);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
    for(int j=0;j<npts_inter;j++){
      // key of gauss point grad value at gauss point
      ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
      _space->gradfuvw(e,uem,vem,0.,Gradm);
      _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      _space->hessfuvw(e,uem,vem,0.,Hessm);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      LocalBasis lbm,lbs;
      lbm.set(e,Gradm,Hessm);
      lbs.set(ie,Grads,lbm.gett0());
      lbm.set_pushForward(&lbs);
      IPState* ips = (*vips)[j];
      IPVariablePlateOnInterface *ipv = dynamic_cast<IPVariablePlateOnInterface*>(ips->getState(ws));
      ipv->setLocalBasis(lbm);
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,nbFF,nbdof,disp,Gradm,Hessm);
      Gradm.clear(); Grads.clear(); Hessm.clear();
    }
    disp.clear();
  }
  // bulk
  for (groupOfElements::elementContainer::const_iterator it = ef->g->begin(); it != ef->g->end(); ++it){
    MElement *e = *it;
    int edge = e->getNumEdges();
    int npts_bulk= ef->getBulkGaussIntegrationRule()->getIntPoints(e,&GP);
    int nbdof = _space->getNumKeys(e);
    int nbFF = e->getNumVertices();
    _ufield->get(e,disp);
    std::vector<IPState*> *vips = _AIPS->getIPstate(e->getNum());
    for(int j=0;j<npts_bulk;j++){
      //grad value at gauss point
      _space->gradfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      _space->hessfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Hessm);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      IPState* ips = (*vips)[j];
      IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
      ipv->setLocalBasis(e,Grads,Hessm);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,nbFF,nbdof,disp,Grads,Hessm);
      // appened method in gradfuvw
      Grads.clear(); Hessm.clear();
    }
    disp.clear();
  }
}

// Compute ipv on a interfaceElement
void IPField::computeIpvShellPlaneStress(MInterfaceElement *ie, const int num, const IPState::whichState ws,
                                          partDomain* ef, const bool virt){
  SVector3 val; // value of a vertex displacement
  IntPt *GP;
  std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
  std::vector<TensorialTraits<double>::HessType> Hessm,Hessp;
  double uem,vem,uep,vep;
  std::vector<double> dispm,dispp;
  dgPartDomain* dgdom = dynamic_cast<dgPartDomain*>(ef);
  if(!virt){
    MElement *em = ie->getElem(0);
    MElement *ep = ie->getElem(1);
    // gauss point
    int npts = dgdom->getInterfaceGaussIntegrationRule()->getIntPoints(ie,&GP);
    // vector with nodal displacement
    int nbdofm = _space->getNumKeys(em);
    int nbdofp = _space->getNumKeys(ep);
    int nbFFm = em->getNumVertices();
    int nbFFp = ep->getNumVertices();
    _ufield->get(em,dispm);
    _ufield->get(ep,dispp);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
    for(int j=num;j<num+npts;j++){
      // key of gauss point grad value at gauss point
      ie->getuvOnElem(GP[j-num].pt[0],uem,vem,uep,vep);
      _space->gradfuvw(em,uem,vem,0.,Gradm);
      _space->hessfuvw(em,uem,vem,0.,Hessm);
      _space->gradfuvw(ep,uep,vep,0.,Gradp);
      _space->hessfuvw(ep,uep,vep,0.,Hessp);
      _space->gradfuvw(ie,GP[j-num].pt[0],GP[j-num].pt[1],GP[j-num].pt[2],Grads);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      LocalBasis lbm,lbp,lbs;
      lbm.set(em,Gradm,Hessm);
      lbp.set(ep,Gradp,Hessp);
      lbs.set(ie,Grads,lbm.gett0(),lbp.gett0());
      lbm.set_pushForward(&lbs);
      lbp.set_pushForward(&lbs);
      IPState* ips = (*vips)[j];
      IPVariablePlateOnInterface *ipv = dynamic_cast<IPVariablePlateOnInterface*>(ips->getState(ws));
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      if(num==0){
        ipv->setLocalBasis(lbm);
        ipv->computeStressAndDeformation(mlaw,nbFFm,nbdofm,dispm,Gradm,Hessm);
      }
      else{
        ipv->setLocalBasis(lbp);
        ipv->computeStressAndDeformation(mlaw,nbFFp,nbdofp,dispp,Gradp,Hessp);
      }

      // appened method in gradfuvw
       Grads.clear(); Gradm.clear();Gradp.clear(); Hessm.clear(); Hessp.clear();
    }
    dispm.clear(); dispp.clear();
  }
  else{  // Virtual interface element
  std::vector<double> disp;
    MElement *e = ie->getElem(0);
    int edge = ie->getEdgeNumber(0);
    int npts_inter= dgdom->getInterfaceGaussIntegrationRule()->getIntPoints(ie,&GP);
    // vector with nodal displacement
    int nbdof = _space->getNumKeys(e);
    int nbFF = e->getNumVertices();
    _ufield->get(e,disp);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
    for(int j=0;j<npts_inter;j++){
      // key of gauss point grad value at gauss point
      ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
      _space->gradfuvw(e,uem,vem,0.,Gradm);
      _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      _space->hessfuvw(e,uem,vem,0.,Hessm);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      LocalBasis lbm,lbs;
      lbm.set(e,Gradm,Hessm);
      lbs.set(ie,Grads,lbm.gett0());
      lbm.set_pushForward(&lbs);
      IPState* ips = (*vips)[j];
      IPVariablePlateOnInterface *ipv = dynamic_cast<IPVariablePlateOnInterface*>(ips->getState(ws));
      ipv->setLocalBasis(lbm);
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,nbFF,nbdof,disp,Gradm,Hessm);
      Gradm.clear(); Grads.clear(); Hessm.clear();
    }
    disp.clear();
  }
}

void IPField::computeIpvShellPlaneStress(MElement *e, IPState::whichState ws, partDomain* ef){
  SVector3 val; // value of a vertex displacement
  IntPt *GP;
  std::vector<TensorialTraits<double>::GradType> Grads;
  std::vector<TensorialTraits<double>::HessType> Hessm;
  std::vector<double> disp;
  // bulk
  int edge = e->getNumEdges();
  int npts_bulk=ef->getBulkGaussIntegrationRule()->getIntPoints(e,&GP);
  int nbdof = _space->getNumKeys(e);
  int nbFF = e->getNumVertices();
  _ufield->get(e,disp);
  std::vector<IPState*> *vips = _AIPS->getIPstate(e->getNum());
  for(int j=0;j<npts_bulk;j++){
    //grad value at gauss point
    _space->gradfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
    _space->hessfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Hessm);
    // local basis on element is needed to compute the local basis on interfaceelement (normal)
    IPState* ips = (*vips)[j];
    IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
    ipv->setLocalBasis(e,Grads,Hessm);
    linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
    ipv->computeStressAndDeformation(mlaw,nbFF,nbdof,disp,Grads,Hessm);
    // appened method in gradfuvw
    Grads.clear(); Hessm.clear();
  }
  disp.clear();
}


// Function to compute VonMises stress
double IPField::getVMShellPlaneStress(MElement *ele, const IPState::whichState ws, const int num,
                                       const partDomain *elas) const{
  double svm =0.;
  if(num<10000){ // VonMises at a Gauss Point
    IPState* ips = (*_AIPS->getIPstate(ele->getNum()))[num];
    IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
    // Stress in plate basis
    double sx,sy,txy,sz,txz,tyz;
    const LocalBasis *lb = ipv->getLocalBasis();
    reductionElement rstress;
    rstress(0,0) = ipv->getSigma(component::xx);
    rstress(1,1) = ipv->getSigma(component::yy);
    rstress(0,1) = rstress(1,0) = ipv->getSigma(component::xy);
    stressTensor stress(rstress,lb);
    // compute of stress for representation
    sz = tyz = txz =0.;
    sx = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
    sy = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
    txy = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
    svm = sqrt(1./2.*((sx-sy)*(sx-sy)+(sy-sz)*(sy-sz)+(sx-sz)*(sx-sz)+6*(txy*txy+tyz*tyz+txz*txz)));
  }
  else{  // Particular value (max,min,mean)
    // loop on IP point
    double svmp;
    IntPt *GP;
    int npts = elas->getBulkGaussIntegrationRule()->getIntPoints(ele,&GP);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
    for(int i=0;i<npts;i++){
      IPState* ips = (*vips)[i];
      IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
      // Stress in plate basis
      double sx,sy,txy,sz,txz,tyz;
      const LocalBasis *lb = ipv->getLocalBasis();
      reductionElement rstress;
      rstress(0,0) = ipv->getSigma(component::xx);
      rstress(1,1) = ipv->getSigma(component::yy);
      rstress(0,1) = rstress(1,0) = ipv->getSigma(component::xy);
      stressTensor stress(rstress,lb);
      // compute of stress for representation
      sz = txz = tyz =0.;
      sx = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
      sy = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
      txy = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
      svmp = sqrt(1./2.*((sx-sy)*(sx-sy)+(sy-sz)*(sy-sz)+(sx-sz)*(sx-sz)+6*(txy*txy+tyz*tyz+txz*txz)));
      if(i==0)
        svm = svmp;
      else{
        switch(num){
          case IPField::max :
            if(svmp>svm) svm=svmp;
            break;
          case IPField::min :
            if(svmp<svm) svm=svmp;
            break;
          case IPField::mean :
            svm+=svmp;
            break;
        }
      }
    }
  if(num==IPField::mean) svm/=(double)npts;
  }
  return svm;
}

// Function to compute VonMises stress
double IPField::getSigmaWithOperationShellPlaneStress(MElement *ele, const IPState::whichState ws, const int num,
                                                      const component::enumcomp cmp, const partDomain *elas) const{
  double sig =0.;
  if(num<10000){ // VonMises at a Gauss Point
    IPState* ips = (*_AIPS->getIPstate(ele->getNum()))[num];
    IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
    // Stress in plate basis
    const LocalBasis *lb = ipv->getLocalBasis();
    reductionElement rstress;
    rstress(0,0) = ipv->getSigma(component::xx);
    rstress(1,1) = ipv->getSigma(component::yy);
    rstress(0,1) = rstress(1,0) = ipv->getSigma(component::xy);
    stressTensor stress(rstress,lb);
    // compute of stress for representation
    switch(cmp){
      case component::xx :
        sig = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
        break;
      case component::yy :
        sig = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
        break;
      case component::xy :
        sig = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
        break;
      default :
        sig=0.;
    }
  }
  else{  // Particular value (max,min,mean)
    // loop on IP point
    double sigp;
    IntPt *GP;
    int npts = elas->getBulkGaussIntegrationRule()->getIntPoints(ele,&GP);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
    for(int i=0;i<npts;i++){
      IPState* ips = (*vips)[i];
      IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
      // Stress in plate basis
      const LocalBasis *lb = ipv->getLocalBasis();
      reductionElement rstress;
      rstress(0,0) = ipv->getSigma(component::xx);
      rstress(1,1) = ipv->getSigma(component::yy);
      rstress(0,1) = rstress(1,0) = ipv->getSigma(component::xy);
      stressTensor stress(rstress,lb);
      switch(cmp){
        case component::xx :
          sigp = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
          break;
        case component::yy :
          sigp = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
          break;
        case component::xy :
          sigp = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
          break;
        default :
          sigp=0.;
      }
      // compute of stress for representation
      if(i==0)
        sig = sigp;
      else{
        switch(num){
          case IPField::max :
            if(sigp>sig) sig=sigp;
            break;
          case IPField::min :
            if(sigp<sig) sig=sigp;
            break;
          case IPField::mean :
            sig+=sigp;
            break;
        }
      }
    }
  if(num==IPField::mean) sig/=(double)npts;
  }
  return sig;
}
