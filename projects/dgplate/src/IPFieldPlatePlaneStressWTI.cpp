//
// C++ Interface: terms
//
// Description: Definition of function defined in class IPField for PlatePlaneStress element with thickness integration
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include"IPField.h"

template<> void IPField<DGelasticField, DgC0FunctionSpace<SVector3> >::compute1statePlatePlaneStressWTI(IPState::whichState ws,
                                                                                                        DGelasticField* ef){
  SVector3 val; // value of a vertex displacement
  IntPt *GP;
  //edge gauss point full dg
  std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
  std::vector<TensorialTraits<double>::HessType> Hessm,Hessp;
  double uem,vem,uep,vep;
  std::vector<double> dispm,dispp;
  for(std::vector<MInterfaceElement*>::iterator it=ef->gi.begin(); it!=ef->gi.end();++it){
    MInterfaceElement *ie = *it;
    MElement *em = ie->getElem(0);
    MElement *ep = ie->getElem(1);
    // gauss point
    int npts = _intBound->getIntPoints(ie,&GP);
    // vector with nodal displacement
    int nbdofm = _space->getNumKeys(em);
    int nbdofp = _space->getNumKeys(ep);
    int nbFFm = em->getNumVertices();
    int nbFFp = ep->getNumVertices();
    dispm.resize(nbdofm);
    _ufield->get(em,dispm);
    dispp.resize(nbdofp);
    _ufield->get(ep,dispp);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
    for(int j=0;j<npts;j++){
      // key of gauss point
      //grad value at gauss point
      ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
      _space->gradfuvw(em,uem,vem,0.,Gradm);
      _space->hessfuvw(em,uem,vem,0.,Hessm);
      _space->gradfuvw(ep,uep,vep,0.,Gradp);
      _space->hessfuvw(ep,uep,vep,0.,Hessp);
      _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      LocalBasis lbm,lbp,lbs;
      lbm.set(em,Gradm);
      lbp.set(ep,Gradp);
      lbs.set(ie,Grads,lbm.gett0(),lbp.gett0());
      lbm.set_pushForward(&lbs);
      lbp.set_pushForward(&lbs);
      IPState* ips = (*vips)[j];
      IPVariablePlateWithThicknessIntegrationOI *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
      ipv->setLocalBasis(lbm);
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,&lbm,nbFFm,nbdofm,dispm,Gradm,Hessm);

      // plus elem
      ips = (*vips)[npts+j];
      ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
      ipv->setLocalBasis(lbp);
      ipv->setLocalBasisOfInterface(lbs);
      ipv->computeStressAndDeformation(mlaw,&lbp,nbFFp,nbdofp,dispp,Gradp,Hessp);

      // appened method in gradfuvw
      Grads.clear(); Gradm.clear();Gradp.clear(); Hessm.clear(); Hessp.clear();
    }
  }

  // Virtual interface element
  std::vector<double> disp;
  for(std::vector<MInterfaceElement*>::iterator it=ef->gib.begin(); it!=ef->gib.end();++it){
    MInterfaceElement *ie = *it;
    // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
    MElement *e = ie->getElem(0);
    int edge = ie->getEdgeNumber(0);
    int npts_inter=_intBound->getIntPoints(ie,&GP);
    // vector with nodal displacement
    int nbdof = _space->getNumKeys(e);
    int nbFF = e->getNumVertices();
    disp.resize(nbdof);
    _ufield->get(e,disp);

    std::vector<IPState*> *vips = _AIPS->getIPstate(ie->getNum());
    for(int j=0;j<npts_inter;j++){
      // key of gauss point
      //grad value at gauss point
      ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
      _space->gradfuvw(e,uem,vem,0.,Gradm);
      _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      _space->hessfuvw(e,uem,vem,0.,Hessm);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      LocalBasis lbm,lbs;
      lbm.set(e,Gradm);
      lbs.set(ie,Grads,lbm.gett0());
      lbm.set_pushForward(&lbs);
      IPState* ips = (*vips)[j];
      IPVariablePlateWithThicknessIntegrationOI *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
      ipv->setLocalBasis(lbm);
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,&lbm,nbFF,nbdof,disp,Gradm,Hessm);
      Gradm.clear(); Grads.clear(); Hessm.clear();
    }
  }
  // bulk
  for (groupOfElements::elementContainer::const_iterator it = ef->g->begin(); it != ef->g->end(); ++it){
    MElement *e = *it;
    int edge = e->getNumEdges();
    int npts_bulk=_intBulk->getIntPoints(e,&GP);
    int nbdof = _space->getNumKeys(e);
    int nbFF = e->getNumVertices();
    disp.resize(nbdof);
    _ufield->get(e,disp);

    std::vector<IPState*> *vips = _AIPS->getIPstate(e->getNum());
    for(int j=0;j<npts_bulk;j++){
      //grad value at gauss point
      _space->gradfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
      _space->hessfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Hessm);
      // local basis on element is needed to compute the local basis on interfaceelement (normal)
      IPState* ips = (*vips)[j];
      IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(ips->getState(ws));
      ipv->setLocalBasis(e,Grads);
      linearElasticLawPlaneStress *mlaw = dynamic_cast<linearElasticLawPlaneStress*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,nbFF,nbdof,disp,Grads,Hessm);
      // appened method in gradfuvw
      Grads.clear(); Hessm.clear();
    }
  }
}

template<> double IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getVMPlatePlaneStressWTI(MElement *ele,
                                                                                                 const IPState::whichState ws,
                                                                                                 const int num,
                                                                                                 const DGelasticField *elas,
                                                                                                 const int pos) const{
  double svm =0.;
  if(num<10000){ // VonMises at a Gauss Point
    IPState* ips = (*_AIPS->getIPstate(ele->getNum()))[num];
    IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(ips->getState(ws));
    // Stress in plate basis
    reductionElement rstress;
    rstress(0,0)=ipv->getSigma(pos,component::xx);
    rstress(1,1)=ipv->getSigma(pos,component::yy);
    rstress(0,1) = rstress(1,0) = ipv->getSigma(pos,component::xy);
    double sx,sy,txy,sz,txz,tyz;
    LocalBasis *lb = ipv->getLocalBasis();
    // get stress tensor from stress reduced
    stressTensor stress(rstress,lb);
    sx = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
    sy = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
    txy = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
    // compute of stress for representation
    sz = tyz = txz =0.;
    svm = sqrt(1./2.*((sx-sy)*(sx-sy)+(sy-sz)*(sy-sz)+(sx-sz)*(sx-sz)+6*(txy*txy+tyz*tyz+txz*txz)));
  }
  else{  // Particular value on element (max,min,mean)
    // loop on IP point
    double svmp;
    IntPt *GP;
    int npts = _intBulk->getIntPoints(ele,&GP);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
	for(int i=0;i<npts;i++){
      IPState* ips = (*vips)[i];
      IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(ips->getState(ws));
      // Stress in plate basis
      double sx,sy,txy,sz,txz,tyz;
      LocalBasis *lb = ipv->getLocalBasis();
      reductionElement rstress;
      rstress(0,0)=ipv->getSigma(pos,component::xx);
      rstress(1,1)=ipv->getSigma(pos,component::yy);
      rstress(0,1) = rstress(1,0) = ipv->getSigma(pos,component::xy);
      // get stress tensor from stress reduced
      stressTensor stress(rstress,lb);
      // compute of stress for representation
      sx = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
      sy = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
      txy = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
      sz = txz = tyz =0.;

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


 template<> const LocalBasis* IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressReducedWTI(MInterfaceElement *iele,
                                                                                                       const IPState::whichState ws,
                                                                                                       const int num,
                                                                                                       const DGelasticField *elas,
                                                                                                       reductionElement &stressTensor,
                                                                                                       const int pos)
{
  if(num<10000){ // stressTensor at a Gauss Point
    IPState* ips = (*_AIPS->getIPstate(iele->getNum()))[num];
    IPVariablePlateWithThicknessIntegrationOI *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
    // Stress in plate basis
    double sa=ipv->getSigma(pos,component::xx), sb=ipv->getSigma(pos,component::yy), sc=0.;
    double tab = ipv->getSigma(pos,component::xy), tac=0. , tbc=0.;
    const LocalBasis *lb = ipv->getLocalBasis();
    // compute of stress for representation
    stressTensor(0,0) = sa;
    stressTensor(0,1) = stressTensor(1,0) = tab;
    stressTensor(1,1) = sb;
    return lb;
  }
  else{  // Particular value on element (max,min,mean)
    // loop on IP point
    stressTensor.setAll(0.);
    IntPt *GP;
    int npts = _intBound->getIntPoints(iele,&GP);
    std::vector<IPState*> *vips = _AIPS->getIPstate(iele->getNum());
    IPVariablePlateWithThicknessIntegrationOI *ipv;
	for(int i=0;i<npts;i++){
      IPState* ips = (*vips)[i];
      ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
      // Stress in plate basis
      double sa=ipv->getSigma(pos,component::xx), sb=ipv->getSigma(pos,component::yy), sc=0.;
      double tab = ipv->getSigma(pos,component::xy), tac=0. , tbc=0.;
      LocalBasis *lb = ipv->getLocalBasis();
      // compute of stress for representation
      if(i==0){
        stressTensor(0,0) = sa;
        stressTensor(0,1) = stressTensor(1,0) = tab;
        stressTensor(1,1) = sb;
      }
      else{
        switch(num){
          case IPField::max :
            Msg::Error("Max on gauss point of stress tensor is not implemented");
            break;
          case IPField::min :
            Msg::Error("Min on gauss point of stress tensor is not implemented");
            break;
          case IPField::mean :
            stressTensor(0,0) += sa;
            stressTensor(0,1) += tab;
            stressTensor(1,1) += sb;
            break;
        }
      }
    }
    if(num==IPField::mean)
      for(int k=0;k<2;k++)
        for(int kk=0;kk<2;kk++) stressTensor(k,kk)/=(double)npts;
    return ipv->getLocalBasis();
  }
}

 template<> const void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressReducedWTI(MInterfaceElement *iele,
                                                                                                       const IPState::whichState ws,
                                                                                                       const int num,
                                                                                                       const DGelasticField *elas,
                                                                                                       reductionElement &stressTensor,
                                                                                                       const int pos, const LocalBasis* lbb[2])
{
  if(num<10000){ // stressTensor at a Gauss Point
    IPState* ips = (*_AIPS->getIPstate(iele->getNum()))[num];
    IPVariablePlateWithThicknessIntegrationOI *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
    // Stress in plate basis
    double sa=ipv->getSigma(pos,component::xx), sb=ipv->getSigma(pos,component::yy), sc=0.;
    double tab = ipv->getSigma(pos,component::xy), tac=0. , tbc=0.;
    const LocalBasis *lb = ipv->getLocalBasis();
    // compute of stress for representation
    stressTensor(0,0) = sa;
    stressTensor(0,1) = stressTensor(1,0) = tab;
    stressTensor(1,1) = sb;
    lbb[0] = lb;
    lbb[1] = ipv->getLocalBasisOfInterface();
  }
  else{  // Particular value on element (max,min,mean)
    // loop on IP point
    stressTensor.setAll(0.);
    IntPt *GP;
    int npts = _intBound->getIntPoints(iele,&GP);
    std::vector<IPState*> *vips = _AIPS->getIPstate(iele->getNum());
    IPVariablePlateWithThicknessIntegrationOI *ipv;
	for(int i=0;i<npts;i++){
      IPState* ips = (*vips)[i];
      ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>(ips->getState(ws));
      // Stress in plate basis
      double sa=ipv->getSigma(pos,component::xx), sb=ipv->getSigma(pos,component::yy), sc=0.;
      double tab = ipv->getSigma(pos,component::xy), tac=0. , tbc=0.;
      LocalBasis *lb = ipv->getLocalBasis();
      // compute of stress for representation
      if(i==0){
        stressTensor(0,0) = sa;
        stressTensor(0,1) = stressTensor(1,0) = tab;
        stressTensor(1,1) = sb;
      }
      else{
        switch(num){
          case IPField::max :
            Msg::Error("Max on gauss point of stress tensor is not implemented");
            break;
          case IPField::min :
            Msg::Error("Min on gauss point of stress tensor is not implemented");
            break;
          case IPField::mean :
            stressTensor(0,0) += sa;
            stressTensor(0,1) += tab;
            stressTensor(1,1) += sb;
            break;
        }
      }
    }
    if(num==IPField::mean)
      for(int k=0;k<2;k++)
        for(int kk=0;kk<2;kk++) stressTensor(k,kk)/=(double)npts;
    lbb[0] = ipv->getLocalBasis();
    lbb[1] = ipv->getLocalBasisOfInterface();
  }
}


template<> double IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressWithOperationPlatePlaneStressWTI(MElement *ele,
                                                                                                 const IPState::whichState ws,
                                                                                                 const int num, const component::enumcomp cmp,
                                                                                                 const DGelasticField *elas,
                                                                                                 const int pos) const{
  double sig =0.;
  if(num<10000){ // VonMises at a Gauss Point
    IPState* ips = (*_AIPS->getIPstate(ele->getNum()))[num];
    IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(ips->getState(ws));
    // Stress in plate basis
    reductionElement rstress;
    rstress(0,0) = ipv->getSigma(pos,component::xx);
    rstress(1,1) = ipv->getSigma(pos,component::yy);
    rstress(0,1) = rstress(1,0) = ipv->getSigma(pos,component::xy);
    LocalBasis *lb = ipv->getLocalBasis();
    stressTensor stress(rstress,lb);
    // compute of stress for representation
    switch(cmp){
      case component::xx :
        sig = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
        break;
      case component::xy :
        sig = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
        break;
      case component::yy :
        sig = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
      default :
        sig=0.;
    }
  }
  else{  // Particular value on element (max,min,mean)
    // loop on IP point
    double sigp;
    IntPt *GP;
    int npts = _intBulk->getIntPoints(ele,&GP);
    std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
	for(int i=0;i<npts;i++){
      IPState* ips = (*vips)[i];
      IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>(ips->getState(ws));
      // Stress in plate basis
      reductionElement rstress;
      rstress(0,0) = ipv->getSigma(pos,component::xx);
      rstress(1,1) = ipv->getSigma(pos,component::yy);
      rstress(0,1) = rstress(1,0) = ipv->getSigma(pos,component::xy);
      LocalBasis *lb = ipv->getLocalBasis();
      stressTensor stress(rstress,lb);
      // compute of stress for representation
      switch(cmp){
        case component::xx :
          sigp = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(0));
          break;
        case component::xy :
          sigp = stress.getComponent(lb->getOrthonormalVector(0),lb->getOrthonormalVector(1));
          break;
        case component::yy :
          sigp = stress.getComponent(lb->getOrthonormalVector(1),lb->getOrthonormalVector(1));
          break;
        default :
          sigp =0.;
      }
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
