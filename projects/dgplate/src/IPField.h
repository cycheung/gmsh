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
template<class T1, class T2>
class IPField{
  protected :
    std::vector<T1>* _efield;
    dofManager<double> *_dm;
    T2* _space;
    QuadratureBase *_intBulk;
    QuadratureBase *_intBound;
    AllIPState *_AIPS;

    // function to compute state depends on element type
    void compute1statePlatePlaneStress(IPState::whichState ws, T1* ef){
      DgC0SolverField<SVector3> SField(_dm, _space); //used to interpolate and evaluate gradient
      SVector3 val; // value of a vertex displacement
      IntPt *GP;
      if(!ef->getFormulation()){ // edge gauss point cg/dg (just minus element)
        std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
        std::vector<TensorialTraits<double>::HessType> Hess;
        double uem,vem,uep,vep;
        fullMatrix<double> disp;
        for(std::vector<MInterfaceElement*>::iterator it=ef->gi.begin(); it!=ef->gi.end();++it){
          MInterfaceElement *ie = *it;
          MElement *e = ie->getElem(0);
          // gauss point
          int npts = _intBound->getIntPoints(ie,&GP);
          // vector with nodal displacement
          int nbdof = _space->getNumKeys(e);
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
            _space->gradfuvw(e,uem,vem,0.,Gradm);
            _space->gradfuvw(ie->getElem(1),uep,vep,0.,Gradp);
            _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
            _space->hessfuvw(e,uem,vem,0.,Hess);
            // local basis on element is needed to compute the local basis on interfaceelement (normal)
            LocalBasis lbm,lbp;
            lbm.set(e,Gradm);
            lbp.set(ie->getElem(1),Gradp);
            IPnum key(e->getNum(),IPnum::createTypeWithTwoInts(ie->getEdgeNumber(0),j));
            IPState* ips = _AIPS->getIPstate(&key);
            IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
            ipv->setLocalBasis(ie,Grads,lbm.gett0(),lbp.gett0());
            ipv->computeStressAndDeformation(ef->getMaterialLaw(),&lbm,nbFF,nbdof,disp,Gradm,Hess); // For interfaceElement it uses
                                                                                                    // the LocalBasis of Element
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
            _space->gradfuvw(em,uem,vem,0.,Gradm);
            _space->hessfuvw(em,uem,vem,0.,Hessm);
            _space->gradfuvw(ep,uep,vep,0.,Gradp);
            _space->hessfuvw(ep,uep,vep,0.,Hessp);
            _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
            // local basis on element is needed to compute the local basis on interfaceelement (normal)
            LocalBasis lbm,lbp;
            lbm.set(em,Gradm);
            lbp.set(ep,Gradp);
            IPnum key(em->getNum(),IPnum::createTypeWithTwoInts(ie->getEdgeNumber(0),j));
            IPState* ips = _AIPS->getIPstate(&key);
            IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
            ipv->setLocalBasis(ie,Grads,lbm.gett0(),lbp.gett0());
            ipv->computeStressAndDeformation(ef->getMaterialLaw(),&lbm,nbFFm,nbdofm,dispm,Gradm,Hessm);

            // plus elem
            IPnum keyp(ep->getNum(),IPnum::createTypeWithTwoInts(ie->getEdgeNumber(1),j));
            ips = _AIPS->getIPstate(&keyp);
            ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
            ipv->setLocalBasis(ie,Grads,lbm.gett0(),lbp.gett0());
            ipv->computeStressAndDeformation(ef->getMaterialLaw(),&lbp,nbFFp,nbdofp,dispp,Gradp,Hessp);

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
      for(std::vector<MInterfaceElement*>::iterator it=ef->gib.begin(); it!=ef->gib.end();++it){
        MInterfaceElement *ie = *it;
        // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
        MElement *e = ie->getElem(0);
        int edge = ie->getEdgeNumber(0);
        int npts_inter=_intBound->getIntPoints(ie,&GP);
        // vector with nodal displacement
        int nbdof = _space->getNumKeys(e);
        int nbFF = e->getNumVertices();
        disp.resize(nbdof,1);
        disp.setAll(0.);
        for(int j=0;j<nbFF;j++){
          SField.getVertexDisplacement(e,val,ef->getFormulation(),j);
          disp(j,0) = val(0);
          disp(j+nbFF,0) = val(1);
          disp(j+2*nbFF,0) = val(2);
        }
        for(int j=0;j<npts_inter;j++){
          IPnum key=IPnum(e->getNum(),IPnum::createTypeWithTwoInts(edge,j));
          // key of gauss point
          //grad value at gauss point
          ie->getuvOnElem(GP[j].pt[0],uem,vem,uep,vep);
          _space->gradfuvw(e,uem,vem,0.,Gradm);
          _space->gradfuvw(ie,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
          _space->hessfuvw(e,uem,vem,0.,Hess);
          // local basis on element is needed to compute the local basis on interfaceelement (normal)
          LocalBasis lbm;
          lbm.set(e,Gradm);
          IPState* ips = _AIPS->getIPstate(&key);
          IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
          ipv->setLocalBasis(ie,Grads,lbm.gett0());
          ipv->computeStressAndDeformation(ef->getMaterialLaw(),&lbm,nbFF,nbdof,disp,Gradm,Hess);
          Gradm.clear(); Grads.clear(); Hess.clear();
        }
      }
      // bulk
      for (groupOfElements::elementContainer::const_iterator it = ef->g->begin(); it != ef->g->end(); ++it){
        MElement *e = *it;
        int edge = e->getNumEdges();
        int npts_bulk=_intBulk->getIntPoints(e,&GP);
        int nbdof = _space->getNumKeys(e);
        int nbFF = e->getNumVertices();
        disp.resize(nbdof,1);
        disp.setAll(0.);
        for(int j=0;j<nbFF;j++){
          SField.getVertexDisplacement(e,val,ef->getFormulation(),j);
          disp(j,0) = val(0);
          disp(j+nbFF,0) = val(1);
          disp(j+2*nbFF,0) = val(2);
        }
        for(int j=0;j<npts_bulk;j++){
          IPnum key=IPnum(e->getNum(),IPnum::createTypeWithTwoInts(edge,j));
          //grad value at gauss point
          _space->gradfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Grads);
          _space->hessfuvw(e,GP[j].pt[0],GP[j].pt[1],GP[j].pt[2],Hess);
          // local basis on element is needed to compute the local basis on interfaceelement (normal)
          IPState* ips = _AIPS->getIPstate(&key);
          IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
          ipv->setLocalBasis(e,Grads);
          ipv->computeStressAndDeformation(ef->getMaterialLaw(),nbFF,nbdof,disp,Grads,Hess);
          // appened method in gradfuvw
          Grads.clear(); Hess.clear();
        }
      }
    }
    double getVMPlatePlaneStress(MElement *ele, const IPState::whichState ws,
                                 const int num, const DGelasticField *elas) const{
      double svm =0.;
      if(num<10000){ // VonMises at a Gauss Point
        IPnum key=IPnum(ele->getNum(),IPnum::createTypeWithTwoInts(ele->getNumEdges(),num)); //modif if on interface (edge instead of getNumVertices)
        IPState* ips = _AIPS->getIPstate(&key);
        IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
        double sx=ipv->getSigma(component::xx), sy=ipv->getSigma(component::yy), sz=0.;
        double txy = ipv->getSigma(component::xy), txz=0. , tyz=0.;
        svm = 1./6.*((sx-sy)*(sx-sy)+(sy-sz)*(sy-sz)+(sx-sz)*(sx-sz)+6*(txy*txy+tyz*tyz+txz*txz));
      }
      else{  // Particular value (max,min,mean)
         // loop on IP point
         double svmp;
         IntPt *GP;
         int npts = _intBulk->getIntPoints(ele,&GP);
         for(int i=0;i<npts;i++){
           IPnum key=IPnum(ele->getNum(),IPnum::createTypeWithTwoInts(ele->getNumEdges(),i)); //modif if on interface (edge instead of getNumVertices)
           IPState* ips = _AIPS->getIPstate(&key);
           IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
           double sx=ipv->getSigma(component::xx), sy=ipv->getSigma(component::yy), sz=0.;
           double txy = ipv->getSigma(component::xy), txz=0. , tyz=0.;
           svmp = 1./6.*((sx-sy)*(sx-sy)+(sy-sz)*(sy-sz)+(sx-sz)*(sx-sz)+6*(txy*txy+tyz*tyz+txz*txz));
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
           if(num==IPField::mean) svm/=(double)npts;
         }
      }
      return svm;
    }
  public :
  enum ElemValue{mean=10001, max=10002, min=10003}; // enum to select a particular value on element
  IPField(std::vector< T1 >* ef,dofManager<double>* pa,T2* sp,
          QuadratureBase *intb, QuadratureBase *intb1, AllIPState* aips) : _efield(ef), _dm(pa), _space(sp),
                                                                           _intBulk(intb), _intBound(intb1), _AIPS(aips){};
  void compute1state(IPState::whichState ws){
    for(int i=0;i<_efield->size();i++){
      switch((*_efield)[i].getSolElemType()){
        case SolElementType::PlatePlaneStress :
          this->compute1statePlatePlaneStress(ws,&(*_efield)[i]);
      }
    }

  }
  // On element only ?? Higher level to pass the associated element (pass edge num)
  double getVonMises(MElement *ele, const IPState::whichState ws, const int num) const{
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
      default :
        Msg::Error("Function getVonMises doesn't exist for element type : %d",ef->getSolElemType());
        svm = 0.;
      }
    return svm;
  }

};

#endif // IPField
