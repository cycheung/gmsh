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

#include "IPState.h"
#include "DgC0PlateSolverField.h"
// IPVariablePlate
void IPVariablePlate::computeStressAndDeformation(linearElasticLawPlaneStress *mtl,const int nbFF, const int nbdof,
                                                  const std::vector<double> &disp,
                                                  const std::vector<TensorialTraits<double>::GradType> &Grads,
                                                  const std::vector<TensorialTraits<double>::HessType> &Hess){
  // Deformation (small deformation plate for now)
  epsilon[0] = epsilongd(0,0,&lb,Grads,disp); // TODO rewrite epsilongd and rhogd more efficiently
  epsilon[1] = epsilongd(1,1,&lb,Grads,disp);
  epsilon[3] = epsilongd(0,1,&lb,Grads,disp);
  epsilon[2] = epsilon[4] = epsilon[5] = 0.;
  rho[0] = rhogd(0,0,&lb,Hess,disp);
  rho[1] = rhogd(1,1,&lb,Hess,disp);
  rho[3] = rhogd(0,1,&lb,Hess,disp);
  rho[2]=rho[4]=rho[5]=0.;
  // stress thanks to material law
  mtl->stress(&lb,epsilon,sigmaMembrane);
  mtl->stress(&lb,rho,sigmaBending);
}

void IPVariablePlate::computeStressAndDeformation(linearElasticLawPlaneStress *mtl,const LocalBasis *lbe,const int nbFF,
                                                  const int nbdof, const std::vector<double> &disp,
                                                  const std::vector<TensorialTraits<double>::GradType> &Grads,
                                                  const std::vector<TensorialTraits<double>::HessType> &Hess){
  // Deformation (small deformation plate for now)
  epsilon[0] = epsilongd(0,0,lbe,Grads,disp); // TODO rewrite epsilongd and rhogd more efficiently
  epsilon[1] = epsilongd(1,1,lbe,Grads,disp);
  epsilon[3] = epsilongd(0,1,lbe,Grads,disp);
  epsilon[2] = epsilon[4] = epsilon[5] = 0.;
  rho[0] = rhogd(0,0,lbe,Hess,disp);
  rho[1] = rhogd(1,1,lbe,Hess,disp);
  rho[3] = rhogd(0,1,lbe,Hess,disp);
  rho[2] = rho[4] = rho[5] =0.;
  // stress thanks to material law
  mtl->stress(lbe,epsilon,sigmaMembrane);
  mtl->stress(lbe,rho,sigmaBending);
}

// IPVariablePlateWithThicknessIntegration
void IPVariablePlateWithThicknessIntegration::computeStressAndDeformation(linearElasticLawPlaneStress *mtl,const int nbFF,
                                                                          const int nbdof, const std::vector<double> &disp,
                                                                          const std::vector<TensorialTraits<double>::GradType> &Grads,
                                                                          const std::vector<TensorialTraits<double>::HessType> &Hess){
  // Deformation (small deformation plate for now)
  // membrane part is the same for all Simpson's point
  double eps11 = epsilongd(0,0,&lb,Grads,disp);
  double eps22 = epsilongd(1,1,&lb,Grads,disp);
  double eps12 = epsilongd(0,1,&lb,Grads,disp);
  // "\xi^3 independant bending part"
  double rho11 = rhogd(0,0,&lb,Hess,disp);
  double rho22 = rhogd(1,1,&lb,Hess,disp);
  double rho12 = rhogd(0,1,&lb,Hess,disp);
  for(int i=0;i<nsimp;i++){
    epsilon[i][0] = eps11 + zsimp[i]*rho11; // TODO rewrite epsilongd and rhogd more efficiently
    epsilon[i][1] = eps22 + zsimp[i]*rho22;
    epsilon[i][3] = eps12 + zsimp[i]*rho12;
    epsilon[i][2] = epsilon[i][4] = epsilon[i][5] = 0.;
    // stress thanks to material law
    mtl->stress(&lb,epsilon[i],sigma[i]);
  }
}

void IPVariablePlateWithThicknessIntegration::computeStressAndDeformation(linearElasticLawPlaneStress *mtl,const LocalBasis *lbe,
                                                                          const int nbFF, const int nbdof,
                                                                          const std::vector<double> &disp,
                                                                          const std::vector<TensorialTraits<double>::GradType> &Grads,
                                                                          const std::vector<TensorialTraits<double>::HessType> &Hess){
  // Deformation (small deformation plate for now)
  // membrane part is the same for all Simpson's point
  double eps11 = epsilongd(0,0,lbe,Grads,disp);
  double eps22 = epsilongd(1,1,lbe,Grads,disp);
  double eps12 = epsilongd(0,1,lbe,Grads,disp);
  // "\xi^3 independant bending part"
  double rho11 = rhogd(0,0,lbe,Hess,disp);
  double rho22 = rhogd(1,1,lbe,Hess,disp);
  double rho12 = rhogd(0,1,lbe,Hess,disp);
  for(int i=0;i<nsimp;i++){
    epsilon[i][0] = eps11 + zsimp[i]*rho11; // TODO rewrite epsilongd and rhogd more efficiently
    epsilon[i][1] = eps22 + zsimp[i]*rho22;
    epsilon[i][3] = eps12 + zsimp[i]*rho12;
    epsilon[i][2] = epsilon[i][4] = epsilon[i][5] = 0.;
    // stress thanks to material law
    mtl->stress(lbe,epsilon[i],sigma[i]);
  }
}
