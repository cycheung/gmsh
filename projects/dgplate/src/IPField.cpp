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
#include "IPField.h"
#include "SimpsonIntegrationRule.h"
#include "materialLaw.h"

// Functions to get IP data
template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStress(const MElement *ele, const int gaussnum,
                                                                                IPState::whichState st, double stress[6]){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>((*vips)[gaussnum]->getState(st));
  ipv->getSigma(stress);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressMembrane(const MElement *ele, const int gaussnum,
                                                                                        IPState::whichState st, double stress[6],
                                                                                        double & h){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>((*vips)[gaussnum]->getState(st));
  ipv->getSigmaMembrane(stress);
  h=ipv->getThickness();
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressBending(const MElement *ele, const int gaussnum,
                                                                                       IPState::whichState st, double stress[6],
                                                                                       double & h){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>((*vips)[gaussnum]->getState(st));
  ipv->getSigmaBending(stress);
  h=ipv->getThickness();
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStress(const MElement *ele, const int gaussnum,
                                                                                IPState::whichState st, double stressM[6],
                                                                                double stressB[6], double & h){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>((*vips)[gaussnum]->getState(st));
  ipv->getSigmaMembrane(stressM);
  ipv->getSigmaBending(stressB);
  h=ipv->getThickness();
}

template<> const LocalBasis* IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressAndLocalBasis(const MElement *ele,
                                                                                                          const int gaussnum,
                                                                                                          IPState::whichState st,
                                                                                                          double stressM[6],
                                                                                                          double stressB[6], double & h){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>((*vips)[gaussnum]->getState(st));
  ipv->getSigmaMembrane(stressM);
  ipv->getSigmaBending(stressB);
  h=ipv->getThickness();
  return ipv->getLocalBasis();
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressAndLocalBasis(const MInterfaceElement *ele,
                                                                                             const int gaussnum,
                                                                                             IPState::whichState st,
                                                                                             double stressM[6], double stressB[6],
                                                                                             double & h, const LocalBasis* lb[3]){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlateOnInterface *ipv = dynamic_cast<IPVariablePlateOnInterface*>((*vips)[gaussnum]->getState(st));
  ipv->getSigmaMembrane(stressM);
  ipv->getSigmaBending(stressB);
  h=ipv->getThickness();
  lb[0] = ipv->getLocalBasis();
  lb[2] = ipv->getLocalBasisOfInterface();
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStress(const MElement *ele, const int gaussnum,
                                                                                IPState::whichState st, std::vector<tab6> &stress){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>((*vips)[gaussnum]->getState(st));
  stress.resize(ipv->getNumSimp());
  ipv->getSigma(stress);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStress(const MElement *ele, const int gaussnum,
                                                                                IPState::whichState st, std::vector<tab6> &stress,
                                                                                std::vector<double> &hsimp){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>((*vips)[gaussnum]->getState(st));
  short int numSimp = ipv->getNumSimp();
  stress.resize(numSimp);
  hsimp.resize(numSimp);
  ipv->getSigma(stress);
  ipv->getSimpsonPoint(hsimp);
}

template<>   LocalBasis* IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressAndLocalBasis(const MElement *ele,
                                                                                                      const int gaussnum,
                                                                                                      IPState::whichState st,
                                                                                                      std::vector<tab6> &stress,
                                                                                                      std::vector<double> &hsimp){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlateWithThicknessIntegration *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegration*>((*vips)[gaussnum]->getState(st));
  short int numSimp = ipv->getNumSimp();
  stress.resize(numSimp);
  hsimp.resize(numSimp);
  ipv->getSigma(stress);
  ipv->getSimpsonPoint(hsimp);
  return ipv->getLocalBasis();
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressAndLocalBasis(const MInterfaceElement *ele,
                                                                                             const int gaussnum,
                                                                                             IPState::whichState st,
                                                                                             std::vector<tab6> &stress,
                                                                                             std::vector<double> &hsimp,
                                                                                             const LocalBasis *lb[3]){
  std::vector<IPState*> *vips = _AIPS->getIPstate(ele->getNum());
  IPVariablePlateWithThicknessIntegrationOI *ipv = dynamic_cast<IPVariablePlateWithThicknessIntegrationOI*>((*vips)[gaussnum]->getState(st));
  short int numSimp = ipv->getNumSimp();
  stress.resize(numSimp);
  hsimp.resize(numSimp);
  ipv->getSigma(stress);
  ipv->getSimpsonPoint(hsimp);
  lb[0] = ipv->getLocalBasis();
  lb[2] = ipv->getLocalBasisOfInterface();
}

// Functions to compute reduction element
template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressReduction(MElement *ele,const int gaussnum,
                                                                                           SolElementType::eltype et,IPState::whichState ws,
                                                                                            reductionElement &nalpha){
    if(et == SolElementType::PlatePlaneStress ){
      double stress[6];
      double h;
      this->getStressMembrane(ele,gaussnum,ws,stress,h);

      // compute nalpha
      nalpha(0,0) = h*stress[0];
      nalpha(0,1) = nalpha(1,0) = h*stress[3];
      nalpha(1,1) = h*stress[1];
    }
    else if(et == SolElementType::PlatePlaneStressWTI ){
      std::vector<tab6> stress ;
      std::vector<double> hsimp;
      this->getStress(ele,gaussnum,ws,stress,hsimp);
      std::vector<double> temp;
      int msize = hsimp.size();
      double h = hsimp[msize-1]-hsimp[0];
      temp.resize(msize);
      for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
      nalpha(0,0) = SimpsonIntegration(temp,h);
      for(int i=0;i<msize;i++) temp[i]=stress[i][3];
      nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
      for(int i=0;i<msize;i++) temp[i]=stress[i][1];
      nalpha(1,1) = SimpsonIntegration(temp,h);
    }
    else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);

  }

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getMomentReduction(MElement *ele,const int gaussnum,
                                                                                           SolElementType::eltype et,IPState::whichState ws,
                                                                                            reductionElement &malpha){
    if(et == SolElementType::PlatePlaneStress ){
      double stress[6];
      double h;
      this->getStressBending(ele,gaussnum,ws,stress,h);

      // compute nalpha
      double hcubdiv12 = h*h*h/12.;
      malpha(0,0) = hcubdiv12*stress[0];
      malpha(0,1) = malpha(1,0) = hcubdiv12*stress[3];
      malpha(1,1) = hcubdiv12*stress[1];
    }
    else if(et == SolElementType::PlatePlaneStressWTI) {
      std::vector<tab6> stress ;
      std::vector<double> zsimp;
      this->getStress(ele,gaussnum,ws,stress,zsimp);
      std::vector<double> temp;
      int msize = zsimp.size();
      temp.resize(msize);
      for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
      malpha(0,0) = SimpsonIntegration(temp,zsimp);
      for(int i=0;i<msize;i++) temp[i]=stress[i][3];
      malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,zsimp);
      for(int i=0;i<msize;i++) temp[i]=stress[i][1];
      malpha(1,1) = SimpsonIntegration(temp,zsimp);
    }
    else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getReduction(MElement *ele,const int gaussnum,
                                                                                           SolElementType::eltype et,IPState::whichState ws,
                                                                                            reductionElement &nalpha, reductionElement &malpha){
    if(et == SolElementType::PlatePlaneStress ){
      double stressMembrane[6];
      double stressBending[6];
      double h;
      this->getStress(ele,gaussnum,ws,stressMembrane,stressBending,h);

      // compute nalpha
      nalpha(0,0) = h*stressMembrane[0];
      nalpha(0,1) = h*stressMembrane[3];
      nalpha(1,1) = h*stressMembrane[1];
      // compute nalpha
      double hcubdiv12 = h*h*h/12.;
      malpha(0,0) = hcubdiv12*stressBending[0];
      malpha(0,1) = malpha(1,0)*hcubdiv12*stressBending[3];
      malpha(1,1) = hcubdiv12*stressBending[1];
    }
    else if(et == SolElementType::PlatePlaneStressWTI){
      std::vector<tab6> stress ;
      std::vector<double> zsimp;
      this->getStress(ele,gaussnum,ws,stress,zsimp);
      std::vector<double> temp;
      int msize = zsimp.size();
      double h = zsimp[msize-1]-zsimp[0];
      temp.resize(msize);
      for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
      nalpha(0,0) = SimpsonIntegration(temp,h);
      malpha(0,0) = SimpsonIntegration(temp,zsimp);
      for(int i=0;i<msize; i++) temp[i]=stress[i][3];
      nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
      malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,zsimp);
      for(int i=0;i<msize; i++) temp[i]=stress[i][1];
      nalpha(1,1) = SimpsonIntegration(temp,h);
      malpha(1,1) = SimpsonIntegration(temp,zsimp);
    }
    else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}

template<> const LocalBasis* IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getReductionAndLocalBasis(MElement *ele,const int gaussnum,
                                                                                           SolElementType::eltype et,IPState::whichState ws,
                                                                                            reductionElement &nalpha, reductionElement &malpha){
    const LocalBasis *lb;
    if(et == SolElementType::PlatePlaneStress ){
      double stressMembrane[6];
      double stressBending[6];
      double h;
      lb = this->getStressAndLocalBasis(ele,gaussnum,ws,stressMembrane,stressBending,h);

      // compute nalpha
      nalpha(0,0) = h*stressMembrane[0];
      nalpha(0,1) = nalpha(1,0) = h*stressMembrane[3];
      nalpha(1,1) = h*stressMembrane[1];
      // compute malpha
      double hcubdiv12 = h*h*h/12.;
      malpha(0,0) = hcubdiv12*stressBending[0];
      malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
      malpha(1,1) = hcubdiv12*stressBending[1];
      return lb;
    }
    else if((et == SolElementType::PlatePlaneStressWTI) or (et == SolElementType::PlatePlaneStressWF)){
      std::vector<tab6> stress ;
      std::vector<double> zsimp;
      lb = this->getStressAndLocalBasis(ele,gaussnum,ws,stress,zsimp);
      std::vector<double> temp;
      int msize = zsimp.size();
      double h = zsimp[msize-1]-zsimp[0];
      temp.resize(msize);
      for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
      nalpha(0,0) = SimpsonIntegration(temp,h);
      malpha(0,0) = SimpsonIntegration(temp,zsimp);
      for(int i=0;i<msize; i++) temp[i]=stress[i][3];
      nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
      malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,zsimp);
      for(int i=0;i<msize; i++) temp[i]=stress[i][1];
      nalpha(1,1) = SimpsonIntegration(temp,h);
      malpha(1,1) = SimpsonIntegration(temp,zsimp);
      return lb;
    }
    else{
      Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
      return NULL;
    }
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getStressReduction(MInterfaceElement *iele, const int gaussnum,
                                                                                         const int numOfGaussPoint,
                                                                                         SolElementType::eltype et, IPState::whichState ws,
                                                                                         reductionElement &nhatmean){
  if(et == SolElementType::PlatePlaneStress){
    double stressMembrane[6];
    double stressBending[6];
    double h;
    const LocalBasis *lb;
    reductionElement nalpha, nhat;

    // minus element
    lb = this->getStressAndLocalBasis(iele,gaussnum,ws,stressMembrane,stressBending,h);
    // nalpha
    nalpha(0,0) = h*stressMembrane[0];
    nalpha(0,1) = nalpha(1,0) = h*stressMembrane[3];
    nalpha(1,1) = h*stressMembrane[1];
    // hat operation and store in nhatmean
    stressReductionHat(nalpha,lb,nhatmean);

    // plus element
    lb = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stressMembrane,stressBending,h);
    // nalpha_p
    nalpha(0,0) = h*stressMembrane[0];
    nalpha(0,1) = nalpha(1,0) = h*stressMembrane[3];
    nalpha(1,1) = h*stressMembrane[1];
    // hat operation
    stressReductionHat(nalpha,lb,nhat);

    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        nhatmean(i,j) += nhat(i,j);
        nhatmean(i,j)*=0.5;
      }
  }
  else if(et == SolElementType::PlatePlaneStressWTI){

    reductionElement nalpha,nhat;
    const LocalBasis *lb;
    std::vector<double> temp;
    std::vector<tab6> stress ;
    std::vector<double> hsimp;

    // minus element
    lb = this->getStressAndLocalBasis(iele,gaussnum,ws,stress,hsimp);
    int msize = hsimp.size();
    double h = hsimp[msize-1]-hsimp[0];
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    nalpha(0,0) = SimpsonIntegration(temp,h);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    nalpha(1,1) = SimpsonIntegration(temp,h);
    // hat operation and store in nhatmean
    stressReductionHat(nalpha,lb,nhatmean);

    // plus element
    lb = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stress,hsimp);
    msize = hsimp.size();
    h = hsimp[msize-1]-hsimp[0];
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    nalpha(0,0) = SimpsonIntegration(temp,h);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    nalpha(1,1) = SimpsonIntegration(temp,h);
    // hat operation
    stressReductionHat(nalpha,lb,nhat);
    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        nhatmean(i,j) += nhat(i,j);
        nhatmean(i,j)*=0.5;
      }

  }
  else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getMomentReduction(MInterfaceElement *iele, const int gaussnum,
                                                                                         const int numOfGaussPoint,
                                                                                         SolElementType::eltype et, IPState::whichState ws,
                                                                                         reductionElement &mhatmean){
  if(et == SolElementType::PlatePlaneStress){
    double stressMembrane[6];
    double stressBending[6];
    double h;
    const LocalBasis *lb;
    reductionElement malpha, mhat;

    // minus element
    lb = this->getStressAndLocalBasis(iele,gaussnum,ws,stressMembrane,stressBending,h);
    // malpha_m
    double hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation and store in mhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // plus element
    lb = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stressMembrane,stressBending,h);
    // nalpha_p
    hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation
    stressReductionHat(malpha,lb,mhat);

    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j) += mhat(i,j);
        mhatmean(i,j)*=0.5;
      }
  }
  else if(et == SolElementType::PlatePlaneStressWTI){

    reductionElement malpha,mhat;
    LocalBasis *lb;
    std::vector<double> temp;
    std::vector<tab6> stress ;
    std::vector<double> hsimp;

    // minus element
    lb = this->getStressAndLocalBasis(iele,gaussnum,ws,stress,hsimp);
    int msize = hsimp.size();
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    malpha(0,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    malpha(1,1) = SimpsonIntegration(temp,hsimp);
    // hat operation and store in nhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // plus element
    lb = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stress,hsimp);
    msize = hsimp.size();
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    malpha(0,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    malpha(1,1) = SimpsonIntegration(temp,hsimp);
    // hat operation
    stressReductionHat(malpha,lb,mhat);
    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j) += mhat(i,j);
        mhatmean(i,j)*=0.5;
      }

  }
  else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}


template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getMomentReductionAndLocalBasis(MInterfaceElement *iele,
                                                                                         const int gaussnum,const int numOfGaussPoint,
                                                                                         SolElementType::eltype et, IPState::whichState ws,
                                                                                         reductionElement &mhatmean,
                                                                                         const LocalBasis* lbb[3]){
  if(et == SolElementType::PlatePlaneStress){
    double stressMembrane[6];
    double stressBending[6];
    double h;
    const LocalBasis *lb;
    reductionElement malpha, mhat;

    // minus element
    this->getStressAndLocalBasis(iele,gaussnum,ws,stressMembrane,stressBending,h,lbb);
    lb = lbb[0];
    // malpha_m
    double hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation and store in mhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // plus element
    lbb[1] = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stressMembrane,stressBending,h);
    lb = lbb[1];
    // nalpha_p
    hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation
    stressReductionHat(malpha,lb,mhat);

    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j) += mhat(i,j);
        mhatmean(i,j)*=0.5;
      }
  }
  else if((et == SolElementType::PlatePlaneStressWTI)){

    reductionElement malpha,mhat;
    const LocalBasis *lb;
    std::vector<double> temp;
    std::vector<tab6> stress ;
    std::vector<double> hsimp;

    // minus element
    this->getStressAndLocalBasis(iele,gaussnum,ws,stress,hsimp,lbb);
    lb = lbb[0];
    int msize = hsimp.size();
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    malpha(0,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    malpha(1,1) = SimpsonIntegration(temp,hsimp);
    // hat operation and store in nhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // plus element
    lbb[1] = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stress,hsimp);
    lb = lbb[1];
    msize = hsimp.size();
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    malpha(0,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    malpha(1,1) = SimpsonIntegration(temp,hsimp);
    // hat operation
    stressReductionHat(malpha,lb,mhat);
    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j) += mhat(i,j);
        mhatmean(i,j)*=0.5;
      }

  }
  else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}


template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getReductionAndLocalBasis(MInterfaceElement *iele,
                                                                                         const int gaussnum,const int numOfGaussPoint,
                                                                                         SolElementType::eltype et, IPState::whichState ws,
                                                                                         reductionElement &nhatmean,
                                                                                         reductionElement &mhatmean,
                                                                                         const LocalBasis* lbb[3]){
  if(et == SolElementType::PlatePlaneStress){
    double stressMembrane[6];
    double stressBending[6];
    double h;
    const LocalBasis *lb;
    reductionElement malpha, mhat;
    reductionElement nalpha, nhat;

    // minus element
    this->getStressAndLocalBasis(iele,gaussnum,ws,stressMembrane,stressBending,h,lbb);
    lb = lbb[0];
    // nalpha_m
    nalpha(0,0) = h*stressMembrane[0];
    nalpha(0,1) = nalpha(1,0) = h*stressMembrane[3];
    nalpha(1,1) = h*stressMembrane[1];
    // hat operation and store in nhatmean
    stressReductionHat(nalpha,lb,nhatmean);
    // malpha_m
    double hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation and store in mhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // plus element
    lbb[1] = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stressMembrane,stressBending,h);
    lb = lbb[1];
    // nalpha_p
    nalpha(0,0) = h*stressMembrane[0];
    nalpha(0,1) = nalpha(1,0) = h*stressMembrane[3];
    nalpha(1,1) = h*stressMembrane[1];
    // hat operation
    stressReductionHat(nalpha,lb,nhat);
    // malpha_p
    hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation
    stressReductionHat(malpha,lb,mhat);

    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j) += mhat(i,j);
        nhatmean(i,j) += nhat(i,j);
        mhatmean(i,j)*=0.5;
        nhatmean(i,j)*=0.5;
      }
  }
  else if((et == SolElementType::PlatePlaneStressWTI) or (et == SolElementType::PlatePlaneStressWF) ){

    reductionElement malpha,mhat;
    reductionElement nalpha,nhat;
    const LocalBasis *lb;
    std::vector<double> temp;
    std::vector<tab6> stress ;
    std::vector<double> hsimp;

    if(!this->getBroken(iele,gaussnum,et,ws)){
      // minus element
      this->getStressAndLocalBasis(iele,gaussnum,ws,stress,hsimp,lbb);
      lb = lbb[0];
      int msize = hsimp.size();
      double h = hsimp[msize-1]-hsimp[0];
      temp.resize(msize);
      for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
      nalpha(0,0) = SimpsonIntegration(temp,h);
      malpha(0,0) = SimpsonIntegration(temp,hsimp);
      for(int i=0;i<msize;i++) temp[i]=stress[i][3];
      nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
      malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
      for(int i=0;i<msize;i++) temp[i]=stress[i][1];
      nalpha(1,1) = SimpsonIntegration(temp,h);
      malpha(1,1) = SimpsonIntegration(temp,hsimp);
      // hat operation and store in nhatmean
      stressReductionHat(nalpha,lb,nhatmean);
      // hat operation and store in nhatmean
      stressReductionHat(malpha,lb,mhatmean);

      // plus element
      lbb[1] = this->getStressAndLocalBasis(iele,gaussnum+numOfGaussPoint,ws,stress,hsimp);
      lb = lbb[1];
      msize = hsimp.size();
      temp.resize(msize);
      for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
      nalpha(0,0) = SimpsonIntegration(temp,h);
      malpha(0,0) = SimpsonIntegration(temp,hsimp);
      for(int i=0;i<msize;i++) temp[i]=stress[i][3];
      nalpha(0,1) = nalpha(1,0) = SimpsonIntegration(temp,h);
      malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
      for(int i=0;i<msize;i++) temp[i]=stress[i][1];
      nalpha(1,1) = SimpsonIntegration(temp,h);
      malpha(1,1) = SimpsonIntegration(temp,hsimp);
      // hat operation
      stressReductionHat(nalpha,lb,nhat);
      stressReductionHat(malpha,lb,mhat);
      // mean
      for(int i=0;i<2;i++)
        for(int j=0;j<2;j++){
          nhatmean(i,j) += nhat(i,j);
          nhatmean(i,j)*=0.5;
          mhatmean(i,j) += mhat(i,j);
          mhatmean(i,j)*=0.5;
        }
    }
    else{
      this->getReductionByCohesiveLawAndLocalBasis(iele,gaussnum,numOfGaussPoint,nhatmean,mhatmean,lbb);
    }
  }
  else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getVirtualMomentReductionAndLocalBasis(MInterfaceElement *iele,
                                                                                         const int gaussnum,const int numOfGaussPoint,
                                                                                         SolElementType::eltype et, IPState::whichState ws,
                                                                                         reductionElement &mhatmean,
                                                                                         const LocalBasis** lbb){
  if(et == SolElementType::PlatePlaneStress){
    double stressMembrane[6];
    double stressBending[6];
    double h;
    const LocalBasis *lb;
    reductionElement malpha;

    // minus element
    this->getStressAndLocalBasis(iele,gaussnum,ws,stressMembrane,stressBending,h,lbb);
    lb = lbb[0];
    // malpha_m
    double hcubdiv12 = h*h*h/12.;
    malpha(0,0) = hcubdiv12*stressBending[0];
    malpha(0,1) = malpha(1,0) = hcubdiv12*stressBending[3];
    malpha(1,1) = hcubdiv12*stressBending[1];
    // hat operation and store in mhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j)*=0.5;
      }
  }
  else if((et == SolElementType::PlatePlaneStressWTI) or (et == SolElementType::PlatePlaneStressWF)){

    reductionElement malpha;
    const LocalBasis *lb;
    std::vector<double> temp;
    std::vector<tab6> stress ;
    std::vector<double> hsimp;

    // minus element
    this->getStressAndLocalBasis(iele,gaussnum,ws,stress,hsimp,lbb);
    lb = lbb[0];
    int msize = hsimp.size();
    temp.resize(msize);
    for(int i=0;i<msize;i++) temp[i]=stress[i][0]; // retrieve component change this ?? TODO
    malpha(0,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][3];
    malpha(0,1) = malpha(1,0) = SimpsonIntegration(temp,hsimp);
    for(int i=0;i<msize;i++) temp[i]=stress[i][1];
    malpha(1,1) = SimpsonIntegration(temp,hsimp);
    // hat operation and store in nhatmean
    stressReductionHat(malpha,lb,mhatmean);

    // mean
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
        mhatmean(i,j)*=0.5;
      }

  }
  else Msg::Error("GetStressReduction is not implemented for SolelemType %d",et);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::getReductionFracture(const MInterfaceElement *iele,
                                                                                             const int numgauss, const int npts,
                                                                                             const SolElementType::eltype elemtype,
                                                                                             const std::vector<double> disp,
                                                                                             const int nbFF_m, const int nbdofm,
                                                                                             const std::vector<TensorialTraits<double>::ValType> &Vals_m,
                                                                                             const std::vector<TensorialTraits<double>::GradType> &Grads_m,
                                                                                             const int nbFF_p,
                                                                                             const std::vector<TensorialTraits<double>::ValType> &Vals_p,
                                                                                             const std::vector<TensorialTraits<double>::GradType> &Grads_p,
                                                                                             reductionElement &nhatmean,
                                                                                             reductionElement &mhatmean) const
{
  // get ipv
  IPVariablePlateOIWF *ipv;
  IPVariablePlateOIWF *ipvp;
  std::vector<IPState*> *vips;
  IPState *ips;
  vips = _AIPS->getIPstate(iele->getNum());
  ips = (*vips)[numgauss];
  ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));
  ips = (*vips)[numgauss+npts];
  ipvp = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::current));

  // Compute the value of delta with the perturbation
  SVector3 ujump;
  double rjump[3];
  displacementjump(Vals_m,nbFF_m,Vals_p,nbFF_p,disp,ujump);
  rotationjump(ipv->getLocalBasisOfInterface(),Grads_m,nbFF_m,nbdofm,ipv->getLocalBasis(),Grads_p,nbFF_p,ipvp->getLocalBasis(),disp,rjump);
  //double delta = ipv->computeDelta(ujump,rjump,ipv->getLocalBasisOfInterface());
  double deltan = ipv->computeDeltaNormal(ujump,rjump,ipv->getLocalBasisOfInterface());
  double deltat = ipv->computeDeltaTangent(ujump,rjump,ipv->getLocalBasisOfInterface());
  // find elasticField (to know the law to use)
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

  // use the law
  linearElasticLawPlaneStressWithFracture *mlaw = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(ef->getMaterialLaw());
  mlaw->getCohesiveReduction(ipv->getm0(),ipv->getn0(),deltan,ipv->getDeltanmax(),deltat,ipv->getDeltatmax(), ipv->getDeltac(),ipv->ifTension(),nhatmean,mhatmean);
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::initialBroken(MInterfaceElement *iele, materialLaw* mlaw ){
  Msg::Info("Interface element %d is broken at initialization",iele->getNum());
  // find the gauss points associated to the interface element iele
  // set broken = true, _tension = true, Gc and sigmac are chosen to the physical value
  // but n0 = m0 = 0 --> no cohesive forces
  reductionElement n0, m0;
  linearElasticLawPlaneStressWithFracture* mlt = dynamic_cast<linearElasticLawPlaneStressWithFracture*>(mlaw);

  // get ipv
  IPVariablePlateOIWF *ipv;
  std::vector<IPState*> *vips;
  IPState *ips;
  vips = _AIPS->getIPstate(iele->getNum());
  double du[3]; du[0] = du[1] = du[2] = 0.;
  double dr[3]; dr[0] = dr[1] = dr[2] = 0.;
  for(int i=0;i<vips->size();i++){
    ips = (*vips)[i];
    ipv = dynamic_cast<IPVariablePlateOIWF*>(ips->getState(IPState::initial));
    ipv->setBroken(mlt->getSigmac(), mlt->getGc(), mlt->getBeta(), n0,m0,du,dr, ipv->getLocalBasisOfInterface(),true);
  }
  ctp.first = iele;
  ctp.second = vips->size()/2;
}

template<> void IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::initialBroken(GModelWithInterface* pModel, std::vector<int> &vnumphys){
  std::vector<MVertex*> vv;
  for(int i=0;i<vnumphys.size();i++){
    // get the vertex associated to the physical entities
    pModel->getMeshVerticesForPhysicalGroup(1,vnumphys[i],vv);
    // find the InterfaceElement associated to these vertex (identify interior node as degree 2 min)
    for(std::vector<DGelasticField>::iterator itfield = _efield->begin(); itfield != _efield->end(); ++itfield){
      for(std::vector<MInterfaceElement*>::iterator it = itfield->gi.begin(); it!=itfield->gi.end(); ++it)
        for(int k=0;k<vv.size();k++)
          if(vv[k] == (*it)->getVertex(2) ) this->initialBroken(*it, itfield->getMaterialLaw());
    }
    vv.clear();
  }
}
