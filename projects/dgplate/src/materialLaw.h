//
// C++ Interface: terms
//
// Description: Define material law
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _MATERIALLAW_H_
#define _MATERIALLAW_H_
// class with all material laws
class materialLaw{
  public :
    enum matname{linearElasticPlaneStress};
    virtual void stress(const double[6],double[6])=0;
};

// class for linear elastic law
class linearElasticLawPlaneStress : public materialLaw{
  protected :
    const double _E; // YoungModulus //Store ??
    const double _nu; // Poisson ratio // Store ??
    const double C11;
    const double C12;
  public :
    linearElasticLawPlaneStress(const double E, const double nu) : _E(E), _nu(nu), C11(E*nu/((1-nu)*(1+nu))), C12(E/(1+nu)){}
    virtual void stress(const double eps[6],double sig[6]){
        sig[0]=(C11+C12)*eps[0]+C12*eps[1];
        sig[1]=C12*eps[0]+(C11+C12)*eps[1];
        sig[2]=0.;
        sig[3]=C12*eps[3];
        sig[4]=0.;
        sig[5]=0.;
    }
};

#endif //_MATERIALLAW_H_
