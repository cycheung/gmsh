//
// C++ Interface: terms
//
// Description: Definition of function defined in class IPField for PlatePlaneStress element with thickness integration and computation of fracture
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include"IPField.h"
#include"IPState.h"

static inline double scaldot(const SVector3 &a, const SVector3 &b){
  double c=0.;
  for(int i=0;i<3;i++)
    c+= a[i]*b[i];
  return c;
}

template<> void IPField<DGelasticField, DgC0FunctionSpace<SVector3> >::compute1statePlatePlaneStressWF(IPState::whichState ws,
                                                                                                        DGelasticField* ef){
  SVector3 val; // value of a vertex displacement
  IntPt *GP;
  //edge gauss point full dg
  std::vector<TensorialTraits<double>::GradType> Grads,Gradm,Gradp;
  std::vector<TensorialTraits<double>::HessType> Hessm,Hessp;
  double uem,vem,uep,vep;
  SVector3 ujump;
  double rjump[3];
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
      IPVariablePlateOIWF *ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(ws));
      IPVariablePlateOIWF *ipvprev = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::previous));
      ipv->setLocalBasis(lbm);
      ipv->setLocalBasisOfInterface(lbs);
      linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
      ipv->computeStressAndDeformation(mlaw,&lbm,nbFFm,nbdofm,dispm,Gradm,Hessm);
      // compute fracture quantities
      if(ipvprev->getBroken()){ // as fracture is compute before nextStep operation
        // compute the jump
        std::vector<TensorialTraits<double>::ValType> Valm;
        _space->fuvw(ie->getElem(0),uem, vem,0., Valm);
        std::vector<TensorialTraits<double>::ValType> Valp;
        _space->fuvw(ie->getElem(1),uep, vep,0., Valp);
        displacementjump(Valm,nbFFm,Valp,nbFFp,dispm,dispp,ujump);
        rotationjump(&lbs,Gradm,nbFFm,&lbm,Gradp,nbFFp,&lbp,dispm,dispp,rjump);
        IPVariablePlateOIWF *ipvprev = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::previous));
        ipv->setFracture(ipvprev,ujump,rjump,&lbs);
        ips = (*vips)[npts+j];
        ipvprev = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::previous));
        ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(ws));
        ipv->setFracture(ipvprev,ujump,rjump,&lbs);
      }

      // plus elem
      ips = (*vips)[npts+j];
      ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(ws));
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

void IPVariablePlateOIWF::setBroken(const double svm, const double Gc, const double beta_, const reductionElement &nalpha,
                                    const reductionElement &malpha, const double du[3], const double dr[3], const LocalBasis *lbs,
                                    const bool tension_)
{
  if(!broken){ // initialize broken at this gauss point
    double hdiv6 = this->getThickness()/6.;
    double hdiva = this->getThickness()/alpha;
    n0 = nalpha;
    m0 = malpha;
    unjump0 = - scaldot(du,lbs->getphi0(1))/lbs->getphi0(1).norm(); // minus because minus normal of LocalBasis of interface = - phi0,1
    unjump =0.;
    utjump0 = scaldot(du,lbs->getphi0d(0))/lbs->getphi0(0).norm();
    rnjump0 = - scaldot(dr,lbs->getphi0(1))/lbs->getphi0(1).norm(); // minus because minus normal of LocalBasis of interface = - phi0,1
    rnjump = 0.;
    rtjump = 0.;
    rtjump0 = scaldot(dr,lbs->getphi0d(0))/lbs->getphi0(0).norm();
    sigmac = svm;   // sigma can be use for tau_max
    deltac = 2*Gc/sigmac;
    beta = beta_;
    deltanmax =0.;
    deltatmax =0.;
    tension = tension_;
    double M0 = malpha(1,1);
    double N0 = nalpha(1,1);
    double N12 = nalpha(0,1);
    double M12 = malpha(0,1);
    double m0abshdiv6, m12abshdiva;
    if(M0<0) {m0abshdiv6 = - M0/hdiv6; rnjump0 = -rnjump0;}
    else m0abshdiv6 = M0/hdiv6;
    etaI = m0abshdiv6/(m0abshdiv6+N0);
    if(M12<0) {m12abshdiva = - M12/hdiv6; rtjump0 = -rtjump0;}
    else m12abshdiva = M12/hdiv6;
    etaII = m12abshdiva/(m12abshdiva+N12);
    deltan0 = etaI*hdiv6*rnjump0 + (1-etaI)*unjump0;
    deltat0 = etaII*hdiva*rtjump0 + (1-etaII)*utjump0;
    broken = true;
  }
}

void IPVariablePlateOIWF::setFracture(const IPVariablePlateOIWF *ipvprev,const double ujump_[3], const double rjump_[3],
                                      const LocalBasis *lbs)
{
  n0 =ipvprev->n0;
  m0 =ipvprev->m0;
  unjump0 = ipvprev->unjump0;
  rnjump0 = ipvprev->rnjump0;
  utjump0 = ipvprev->utjump0;
  rtjump0 = ipvprev->rtjump0;
  sigmac = ipvprev->sigmac;
  deltac = ipvprev->deltac;
  deltan0 = ipvprev->deltan0;
  deltat0 = ipvprev->deltat0;
  deltanmax = ipvprev->deltanmax;
  deltatmax = ipvprev->deltatmax;
  etaI = ipvprev->etaI;
  etaII= ipvprev->etaII;
  alpha = ipvprev->alpha;
  beta = ipvprev->beta;
  broken = ipvprev->broken;
  tension = ipvprev->tension;
  unjump = - scaldot(ujump_,lbs->getphi0(1))/lbs->getphi0(1).norm() - unjump0; // minus because normal of LocalBasis of interface = - phi0,1
  utjump = scaldot(ujump_,lbs->getphi0(0))/lbs->getphi0(0).norm() - utjump0;
  if( m0(1,1)<0) rnjump = scaldot(rjump_,lbs->getphi0(1))/lbs->getphi0(1).norm() - rnjump0;
  else rnjump = - scaldot(rjump_,lbs->getphi0(1))/lbs->getphi0(1).norm() - rnjump0;
  if( m0(0,1)<0) rtjump = scaldot(rjump_,lbs->getphi0(0))/lbs->getphi0(0).norm() - rtjump0;
  else rtjump = - scaldot(rjump_,lbs->getphi0(0))/lbs->getphi0(0).norm() - rtjump0;
  double hdiv6 = this->getThickness()/6.;
  double hdiva = this->getThickness()/alpha;
  deltan = etaI*hdiv6*rnjump + (1-etaI)*unjump;
  deltat = etaII*hdiva*rtjump+(1-etaII)*utjump;
}

double IPVariablePlateOIWF::computeDelta(const SVector3 ujump_, const double rjump_[3], const LocalBasis *lbs) const
{
  double unor = - scaldot(ujump_,lbs->getphi0(1))/lbs->getphi0(1).norm(); // minus because normal of LocalBasis of interface = - phi0,1
  double ut = scaldot(ujump_,lbs->getphi0d(0))/lbs->getphi0(0).norm();
  double rnor =  - scaldot(rjump_,lbs->getphi0(1))/lbs->getphi0(1).norm();
  double rt = scaldot(rjump_,lbs->getphi0d(0))/lbs->getphi0(0).norm();
  double deltan, deltat;
  if(m0(1,1)>0) deltan = (1-etaI)*(unor-unjump0) + etaI*this->getThickness()/6.*(-rnor-rnjump0);
  else deltan = (1-etaI)*(unor-unjump0) + etaI*this->getThickness()/6.*(rnor-rnjump0);
  if(m0(0,1)>0) deltat = (1-etaII)*(ut-utjump0) + etaII*this->getThickness()/alpha*(-rt-rtjump0);
  else deltat = (1-etaII)*(ut-utjump0) + etaII*this->getThickness()/alpha*(rt-rtjump0);
  return sqrt(deltan*deltan+beta*beta*deltat*deltat);
}

double IPVariablePlateOIWF::computeDeltaNormal(const SVector3 ujump_, const double rjump_[3], const LocalBasis *lbs) const
{
  double unor = - scaldot(ujump_,lbs->getphi0(1))/lbs->getphi0(1).norm(); // minus because normal of LocalBasis of interface = - phi0,1
  double rnor = - scaldot(rjump_,lbs->getphi0(1))/lbs->getphi0(1).norm();
  if(m0(1,1)<0) return (1-etaI)*(unor-unjump0) + etaI*this->getThickness()/6.*(-rnor-rnjump0);
  else return (1-etaI)*(unor-unjump0) + etaI*this->getThickness()/6.*(rnor-rnjump0);
}
double IPVariablePlateOIWF::computeDeltaTangent(const SVector3 ujump_, const double rjump_[3], const LocalBasis *lbs) const
{
  double utan = - scaldot(ujump_,lbs->getphi0(0))/lbs->getphi0(0).norm(); // minus because normal of LocalBasis of interface = - phi0,1
  double rtan = - scaldot(rjump_,lbs->getphi0(0))/lbs->getphi0(0).norm();
  if(m0(0,1)<0) return (1-etaII)*(utan-utjump0) + etaII*this->getThickness()/alpha*(-rtan-rtjump0);
  else return (1-etaII)*(utan-utjump0) + etaII*this->getThickness()/alpha*(rtan-rtjump0);
}
