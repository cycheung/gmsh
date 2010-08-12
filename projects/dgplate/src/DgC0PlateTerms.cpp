//
// C++ Interface: terms
//
// Description: Elementary matrix terms for C0 Dg Plate
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "DgC0PlateTerms.h"
#include "DgC0PlateElementaryTerms.h" // contains some function used in this file to compute elementary matrix
void IsotropicElasticStiffBulkTermC0Plate::get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m)
{
  if (sym)
  {
    // Initialization of some data
    const int nbdof = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(ele);
    const int nbFF = nbdof/3;
    double Cmt= 0., Cnt=0.;
    LinearElasticShellHookeTensor HOOKe; LinearElasticShellHookeTensor *H=&HOOKe; // make a pointer in an other way
    m.resize(nbdof, nbdof);
    m.setAll(0.);
    std::vector<TensorialTraits<double>::HessType> Hess;
    std::vector<TensorialTraits<double>::GradType> Grad;
    double Bn[256][3][2][2], Bm[256][3][2][2]; // max order 256 or dynamical allocation ?? better than std::vector and fullMatrix ?? create a class with this variable
    LocalBasis LB;
    LocalBasis *lb=&LB; // two last line in 1 operation ??
    double me_bending[3][3];
    double me_membrane[3][3]; // better than fullMatrix<double> used for me_bending ??
    // sum on Gauss' points
    for (int i = 0; i < npts; i++)
    {
      // Coordonate of Gauss' point i
      const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
      // Weight of Gauss' point i and Jacobian's value at this point
      const double weight = GP[i].weight;
      // Compute of Hessian of SF. Each component of the vector in link to a shape function.
      // It give the value of second derivative of SF in the isoparametric configuration
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(ele,u, v, w, Grad); // a optimiser the jacobian cannot be given in argument...
      DgC0BilinearTerm<SVector3,SVector3>::space1.hessfuvw(ele,u, v, w, Hess); // a optimiser

      lb->set(ele,Grad,Hess); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> qui retourne un vecteur vide prob de template??

      // multiplication of constant by detJ and weight
      Cmt = lb->getJacobian() * weight *Cm;
      Cnt = lb->getJacobian() * weight *Cn;
      // compute of Hooke tensor
      H->set(lb,nu);

      // compute Bn value
      compute_Bn(lb,Grad,nbFF,Bn);
      Compute_Bm(lb,Grad,Hess,nbFF,Bm);

     // loop on SF to construct the elementary (bulk) stiffness matrix at the Gauss' point i
     for(int j=0; j<nbFF;j++)
       for(int k=0;k<nbFF;k++){
         BulkC0PlateDGStiffnessBendingTerms(Bm[j],Bm[k],H,me_bending);
         BulkC0PlateDGStiffnessMembraneTerms(Bn[j],Bn[k],H,me_membrane);
         for(int jj=0;jj<3;jj++)
           for(int kk=0;kk<3;kk++)
             m(j+(jj*nbFF),k+(kk*nbFF)) += (Cmt*me_bending[jj][kk]+Cnt*me_membrane[jj][kk]);
       }
     // clear the hessian and Grad because the components append in hessfuvw and gradfuvw
     Hess.clear(); Grad.clear();
    }
/*   m.print("bulk");
   m.setAll(0.);
    // By numerical perturbation (Verification OK)
    double eps=1.e-6;
    double epsm = -eps;
    fullVector<double> fp(nbdof);
    fullVector<double> fm(nbdof);
    fp.scale(0.);
    fm.scale(0.);
    //fullMatrix<double> m2(nbdof,nbdof); m2.setAll(0.);
    DgC0LinearTerm<SVector3> *lterm = this->getLinearTerm();
    Dof D(0,0);
    //displacementField ufield(pAssembler,elasticFields,3,LagSpace->getId())
    for(int j=0;j<ele->getNumVertices();j++){
      for(int jj=0;jj<3;jj++){
        if(!fullDg)
          D = Dof(ele->getVertex(j)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000));
        else
          D = Dof(ele->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000,j));
        fm.scale(0.);fp.scale(0.);
        ufield->set(D,eps);
        //ipf->compute1state(IPState::current);
        ipf->computeIpv(ele,IPState::current);
        lterm->get(ele,npts,GP,fp);
        ufield->set(D,epsm);
        ufield->set(D,epsm);
//        ipf->compute1state(IPState::current);
        ipf->computeIpv(ele,IPState::current);
        lterm->get(ele,npts,GP,fm);
        ufield->set(D,eps);
        for(int k=0;k<nbdof;k++)
          m(k,jj*ele->getNumVertices()+j)=(fp(k)-fm(k))/(2.*eps);
      }
    }
//    ipf->compute1state(IPState::current);
      ipf->computeIpv(ele,IPState::current);
    m.print("Bulk pert");*/

  }
  else
    printf("not implemented\n");
}

void IsotropicElasticForceBulkTermC0Plate::get(MElement *ele,int npts,IntPt *GP, fullVector<double> &m)
{
  // Initialization of some data
  const int nbdof = DgC0LinearTerm<SVector3>::space1.getNumKeys(ele);
  const int nbFF = nbdof/3;
  //double Cmt= 0.,Cnt=0.;
  double wJ;
  LinearElasticShellHookeTensor HOOKe; LinearElasticShellHookeTensor *H=&HOOKe; // make a pointer in an other way
  m.resize(nbdof);
  m.scale(0.);
  std::vector<TensorialTraits<double>::HessType> Hess;
  std::vector<TensorialTraits<double>::GradType> Grad;
  double Bm[256][3][2][2];
  reductionElement nalpha,malpha;

  // sum on Gauss' points
  for (int i = 0; i < npts; i++)
  {
    // Coordonate of Gauss' point i
    const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
    // Weight of Gauss' point i and Jacobian's value at this point
    const double weight = GP[i].weight;

    // Compute of Hessian of SF. Each component of the vector in link to a shape function.
    // It give the value of second derivative of SF in the isoparametric configuration
    DgC0LinearTerm<SVector3>::space1.gradfuvw(ele,u, v, w, Grad); // a optimiser the jacobian cannot be given in argument...
    DgC0LinearTerm<SVector3>::space1.hessfuvw(ele,u, v, w, Hess); // a optimiser
    // not very elegant but very usefull
    const LocalBasis *lb = ipf->getReductionAndLocalBasis(ele,i,_elemtype,IPState::current,nalpha,malpha);

    Compute_Bm(lb,Grad,Hess,nbFF,Bm);
    wJ = lb->getJacobian() * weight;
    for(int j=0; j<nbFF;j++){
      for(int k=0;k<3;k++)
        for(int alpha=0;alpha<2;alpha++)
          for(int beta=0;beta<2;beta++)
            m(j+k*nbFF)+=wJ*(Grad[j](alpha)*nalpha(alpha,beta)*lb->getphi0(beta,k)+ Bm[j][k][alpha][beta]*malpha(alpha,beta));
    }
  // clear the hessian and Grad because the components append in hessfuvw and gradfuvw
  Hess.clear(); Grad.clear();
  }
  if(inverseSign)
    for(int i=0;i<nbdof;i++) m(i) = -m(i);
}

void IsotropicElasticStiffInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP, fullMatrix<double> &m)
{
  if (sym)
  {
    // data initialization
    // Retrieve of the element link to the interface element velem[0] == minus elem ; velem == plus elem
    MInterfaceElement *iele = dynamic_cast<MInterfaceElement*>(ele);
    MElement ** velem = iele->getElem();
    const int nbdof_m = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(velem[0]);
    const int nbdof_p = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(velem[1]);
    const int nbFF_m = nbdof_m/3;
    const int nbFF_p = nbdof_p/3;
    // Initialization
    m.resize(nbdof_m+nbdof_p, nbdof_m+nbdof_p);
    m.setAll(0.);
    double uem,uep,vem,vep;
    std::vector<TensorialTraits<double>::ValType> Val_m, Val_p;
    std::vector<TensorialTraits<double>::GradType> Grads_m;
    std::vector<TensorialTraits<double>::GradType> Grads_p;
    std::vector<TensorialTraits<double>::GradType> Grads;
    std::vector<TensorialTraits<double>::HessType> Hess_m;
    std::vector<TensorialTraits<double>::HessType> Hess_p;
    LocalBasis LBS,LBP,LBM;
    LocalBasis *lbs=&LBS; LocalBasis *lb_p=&LBP; LocalBasis *lb_m=&LBM;
    double Bhat_p[256][3][2][2], Bhat_m[256][3][2][2], Bn_m[256][3][2][2], Bn_p[256][3][2][2],Bm_m[256][3][2][2], Bm_p[256][3][2][2];
    double Bs_m[256][3], Bs_p[256][3];
    LinearElasticShellHookeTensor HOOKEhat_p;
    LinearElasticShellHookeTensor HOOKehat_m;
    LinearElasticShellHookeTensor HOOKehatmean;
    LinearElasticShellHookeTensor *Hhat_p=&HOOKEhat_p;
    LinearElasticShellHookeTensor *Hhat_m=&HOOKehat_m;
    LinearElasticShellHookeTensor *Hhatmean=&HOOKehatmean;
    double Deltat_m[256][3][3], Deltat_p[256][3][3];
    double Deltau_m[256][2][3], Deltau_p[256][2][3];
    double Cmt,Cnt,Cst;
    double me_cons[3][3];
    double me_comp[3][3];
    double me_stab[3][3];
    std::vector<bool> vbroken;
    std::vector<bool> vDeltanNeg;
    ipf->getBroken(iele,_elemtype,vbroken,vDeltanNeg);
    // Characteristic size of interface element
    double h_s = iele->getCharacteristicSize();
    const double Bhs = beta1/h_s;
    const double B2hs= beta2/h_s;
    const double B3hs= beta3/h_s;
    // sum on Gauss' points
    for (int i = 0; i < npts; i++)
    {
      // Coordonate of Gauss' point i
      const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
      // Weight of Gauss' point i
      const double weight = GP[i].weight;
      //printf("Abscisse of gauss point %f\n",u);
      //Compute the position (u,v) in the element of the Gauss point (to know where evaluate the shape functions)
      iele->getuvOnElem(u,uem,vem,uep,vep);
      // Compute of gradient and hessian of shape functions on interface element and on elements
      // ATTENTION after multiplication by multipliers (functionspace 276) components are in the "second line of tensor change this ??
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(iele,u, v, w, Grads); // grad of shape fonction on interface element
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(velem[0],uem, vem, w, Grads_m); // w = 0
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(velem[1],uep, vep, w, Grads_p);
      DgC0BilinearTerm<SVector3,SVector3>::space1.hessfuvw(velem[0],uem, vem, w, Hess_m);
      DgC0BilinearTerm<SVector3,SVector3>::space1.hessfuvw(velem[1],uep, vep, w, Hess_p);

      // basis of elements and interface element
      lb_m->set(velem[0],Grads_m,Hess_m); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lb_p->set(velem[1],Grads_p,Hess_p); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lbs->set(iele,Grads,lb_p->gett0(),lb_m->gett0()); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      // PushForwardTensor
      lb_m->set_pushForward(lbs);
      lb_p->set_pushForward(lbs);

      if(!vbroken[i]){
        // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
        Compute_Bm(lb_p,Grads_p,Hess_p,nbFF_p,Bm_p);
        Compute_Bm(lb_m,Grads_m,Hess_m,nbFF_m,Bm_m);
        compute_Bhat(lb_p,nbFF_p,Bm_p,Bhat_p);
        compute_Bhat(lb_m,nbFF_m,Bm_m,Bhat_m);

        // Compute of Hooke hat tensor on minus and plus element
        Cmt = weight * lbs->getJacobian()  * Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
        Hhat_p->hat(lb_p,Cmt,nu);
        Hhat_m->hat(lb_m,Cmt,nu); //Hhat_m->hat(lb_m,H_m);
        Hhatmean->mean(Hhat_p,Hhat_m); // mean value of tensor by component used for stability term

        // Compute of Deltat_tilde
        compute_Deltat_tilde(lb_p,Grads_p,nbFF_p,Deltat_p);
        compute_Deltat_tilde(lb_m,Grads_m,nbFF_m,Deltat_m);
        for(int j=0;j<nbFF_m;j++){
          for(int k=0;k<nbFF_m;k++){
            consC0PlateStiffnessTerms(Hhat_m,Bhat_m[k],Deltat_m[j],lbs,me_cons);
            compC0PlateStiffnessTerms(Hhat_m,Bhat_m[j],Deltat_m[k],lbs,me_comp);
            stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_m[j], Deltat_m[k],lbs,me_stab);
            for(int jj=0;jj<3;jj++)
              for(int kk=0;kk<3;kk++)
                m(j+jj*nbFF_m,k+kk*nbFF_m) += (- me_cons[jj][kk] - me_comp[jj][kk] + Bhs * me_stab[jj][kk] );
          }
          for(int k=0;k<nbFF_p;k++){
            consC0PlateStiffnessTerms(Hhat_p,Bhat_p[k],Deltat_m[j],lbs,me_cons);
            compC0PlateStiffnessTerms(Hhat_m,Bhat_m[j],Deltat_p[k],lbs,me_comp);
            stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_m[j], Deltat_p[k],lbs,me_stab);
            for(int jj=0;jj<3;jj++)
              for(int kk=0;kk<3;kk++)
                m(j+jj*nbFF_m,k+nbdof_m+kk*nbFF_p) += ( - me_cons[jj][kk] + me_comp[jj][kk]- Bhs * me_stab[jj][kk] );
          }
        }
        for(int j=0;j<nbFF_p;j++){
          for(int k=0;k<nbFF_m;k++){
            consC0PlateStiffnessTerms(Hhat_m,Bhat_m[k],Deltat_p[j],lbs,me_cons);
            compC0PlateStiffnessTerms(Hhat_p,Bhat_p[j],Deltat_m[k],lbs,me_comp);
            stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_p[j], Deltat_m[k],lbs,me_stab);
            for(int jj=0;jj<3;jj++)
              for(int kk=0;kk<3;kk++)
                m(j+(jj*nbFF_p)+nbdof_m,k+kk*nbFF_m) += (me_cons[jj][kk] - me_comp[jj][kk] - Bhs * me_stab[jj][kk]);
          }
          for(int k=0;k<nbFF_p;k++){
            consC0PlateStiffnessTerms(Hhat_p,Bhat_p[k],Deltat_p[j],lbs,me_cons);
            compC0PlateStiffnessTerms(Hhat_p,Bhat_p[j],Deltat_p[k],lbs,me_comp);
            stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_p[j], Deltat_p[k],lbs,me_stab);
            for(int jj=0;jj<3;jj++)
              for(int kk=0;kk<3;kk++)
                m(j+(jj*nbFF_p)+nbdof_m,k+(kk*nbFF_p)+nbdof_m) += ( me_cons[jj][kk] + me_comp[jj][kk] + Bhs * me_stab[jj][kk]);
          }
        }
        // If Full Dg formulation extra terms must be added
        if(fullDg){
          // Shape functions evaluated in u,v
          DgC0BilinearTerm<SVector3,SVector3>::space1.fuvw(velem[0],uem, vem, w, Val_m);
          DgC0BilinearTerm<SVector3,SVector3>::space1.fuvw(velem[1],uep, vep, w, Val_p);

          // Compute Bnhat vector (used previous Bhat_m and Bhat_p)
          // first value of Bn are needed
          compute_Bn(lb_m,Grads_m,nbFF_m,Bn_m);
          compute_Bn(lb_p,Grads_p,nbFF_p,Bn_p);
          compute_Bhat(lb_m,nbFF_m,Bn_m,Bhat_m);
          compute_Bhat(lb_p,nbFF_p,Bn_p,Bhat_p);

          // Compute Hooke tensor
          Cnt = weight*lbs->getJacobian()*Cn;
          Hhat_m->hat(lb_m,Cnt,nu); // Redondant avec + haut (faut enlever le Cmt pour le calcul augmente cout calcul ??)
          Hhat_p->hat(lb_p,Cnt,nu);
          Hhatmean->mean(Hhat_m,Hhat_p);

          // Assembly (regroup all ??)
          for(int j=0;j<nbFF_m;j++){
            for(int k=0;k<nbFF_m;k++){
              consC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_m[k],Val_m[j],lbs,me_cons);
              compC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_m[j],Val_m[k],lbs,me_comp);
              stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_m[j], Val_m[k],lbs,me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  m(j+jj*nbFF_m,k+kk*nbFF_m) += (- me_cons[jj][kk] - me_comp[jj][kk] + B2hs * me_stab[jj][kk] );
            }
            for(int k=0;k<nbFF_p;k++){
              consC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_p[k],Val_m[j],lbs,me_cons);
              compC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_m[j],Val_p[k],lbs,me_comp);
              stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_m[j], Val_p[k],lbs,me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  m(j+jj*nbFF_m,k+nbdof_m+kk*nbFF_p) += ( - me_cons[jj][kk] + me_comp[jj][kk]- B2hs * me_stab[jj][kk] );
            }
          }
          for(int j=0;j<nbFF_p;j++){
            for(int k=0;k<nbFF_m;k++){
              consC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_m[k],Val_p[j],lbs,me_cons);
              compC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_p[j],Val_m[k],lbs,me_comp);
              stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_p[j], Val_m[k],lbs,me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  m(j+(jj*nbFF_p)+nbdof_m,k+kk*nbFF_m) += (me_cons[jj][kk] - me_comp[jj][kk] - B2hs * me_stab[jj][kk]);
            }
            for(int k=0;k<nbFF_p;k++){
              consC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_p[k],Val_p[j],lbs,me_cons);
              compC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_p[j],Val_p[k],lbs,me_comp);
              stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_p[j], Val_p[k],lbs,me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  m(j+(jj*nbFF_p)+nbdof_m,k+(kk*nbFF_p)+nbdof_m) += ( me_cons[jj][kk] + me_comp[jj][kk] + B2hs * me_stab[jj][kk]);
            }
          }
          // Out of Plane term
          Cst = weight*lbs->getJacobian()*Cs*B3hs; // Hooke tensor in "shearing" 1 component = Cs

          // compute B vector
          compute_Bs(lbs,Val_m,nbFF_m,Bs_m);
          compute_Bs(lbs,Val_p,nbFF_p,Bs_p);
          for(int j=0;j<nbFF_m;j++){
            for(int k=0;k<nbFF_m;k++){
              stabilityC0PlateStiffnessShearingTerms(Bs_m[j], Bs_m[k],me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  {m(j+jj*nbFF_m,k+kk*nbFF_m) += Cst * me_stab[jj][kk]; }//printf("-- %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);
            }
            for(int k=0;k<nbFF_p;k++){
              stabilityC0PlateStiffnessShearingTerms(Bs_m[j], Bs_p[k],me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  {m(j+jj*nbFF_m,k+nbdof_m+kk*nbFF_p) += -Cst * me_stab[jj][kk]; }//printf("-+ %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);
            }
          }
          for(int j=0;j<nbFF_p;j++){
            for(int k=0;k<nbFF_m;k++){
              stabilityC0PlateStiffnessShearingTerms(Bs_p[j], Bs_m[k],me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  {m(j+(jj*nbFF_p)+nbdof_m,k+kk*nbFF_m) += -Cst * me_stab[jj][kk]; }//printf("+- %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);
            }
            for(int k=0;k<nbFF_p;k++){
              stabilityC0PlateStiffnessShearingTerms(Bs_p[j], Bs_p[k],me_stab);
              for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                  {m(j+(jj*nbFF_p)+nbdof_m,k+(kk*nbFF_p)+nbdof_m) += Cst * me_stab[jj][kk]; }//printf("++ %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);
            }
          }
          Val_m.clear();Val_p.clear();
        }
      }
      else{
        // get matrix by perturbation
        // WARNING : as the distinction between cg/dg and full dg formulation is only made for the assembling
        // the perturbation must be applied on cg/dg case as it is a full dg formulation
        DgC0BilinearTerm<SVector3,SVector3>::space1.fuvw(velem[0],uem, vem, w, Val_m);
        DgC0BilinearTerm<SVector3,SVector3>::space1.fuvw(velem[1],uep, vep, w, Val_p);
        //init
        fullVector<double> fp;
        fullVector<double> fm;
        fp.resize(nbdof_m+nbdof_p);
        fm.resize(nbdof_m+nbdof_p);

        // get displacement vector
        std::vector<double> disp;
        ufield->get(iele,disp);

        // pertubation on displacement vector
        for(int j=0;j<disp.size();j++){
          fm.scale(0.);
          disp[j] -= perturbation;
          // compute fm with a function for now FIX IT TODO ??
          computeFintFrac(iele,disp,lb_m,lb_p,lbs, weight,i,npts,_elemtype,ipf,nbFF_m,Val_m,Grads_m,nbdof_m,nbFF_p,Val_p,Grads_p,fm);
          disp[j]+=perturbation;
          fp.scale(0.);
          disp[j]+=perturbation;
          // compute fp with a function for now FIX IT TODO ??
          computeFintFrac(iele,disp,lb_m,lb_p,lbs, weight,i,npts,_elemtype,ipf,nbFF_m,Val_m,Grads_m,nbdof_m,nbFF_p,Val_p,Grads_p,fp);
          disp[j]-=perturbation;
          // add in the matrix
          for(int k=0;k<nbdof_m+nbdof_p;k++)
            m(k,j) += (fp(k)-fm(k))/(perturbation+perturbation);
        }
        Val_m.clear(); Val_p.clear();
      }
      // Because component are push_back in Grads in gradfuvw idem for hess
      Grads_m.clear(); Grads_p.clear(); Hess_m.clear(); Hess_p.clear(); Grads.clear();
    }
//    m.print("Interface");
/*    m.setAll(0.);
    // By numerical perturbation (Verification OK)
    double eps=1.e-6;
    double per;
    fullVector<double> fp(nbdof_m+nbdof_p);
    fullVector<double> fm(nbdof_m+nbdof_p);
    fp.scale(0.);
    fm.scale(0.);
    //fullMatrix<double> m2(nbdof_m+nbdof_p,nbdof_m+nbdof_p); m2.setAll(0.);
    IsotropicElasticForceInterfaceTermC0Plate *lterm = this->getLinearTerm();
    lterm->setMinus(true);
    Dof D(0,0);
    for(int j=0;j<iele->getElem(0)->getNumVertices();j++){
      for(int jj=0;jj<3;jj++){
        if(!fullDg)
          D = Dof(iele->getElem(0)->getVertex(j)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000));
        else
          D = Dof(iele->getElem(0)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000,j));
        fm.scale(0.);fp.scale(0.);
        per = eps;
        lterm->setPertDof(D);
        lterm->setPert(per);
        ufield->set(D,per);
        ipf->computeIpv(iele,0,IPState::current); // 0 for - elem and npts for + elem
        lterm->get(iele,npts,GP,fp);
        per = -2.*eps;
        ufield->set(D,per);
        lterm->setPert(-eps);
        ipf->computeIpv(iele,0,IPState::current); // 0 for - elem and npts for + elem
        lterm->get(iele,npts,GP,fm);
        per = eps;
        ufield->set(D,per);
        for(int k=0;k<nbdof_m+nbdof_p;k++)
          m(k,jj*iele->getElem(0)->getNumVertices()+j)=(fp(k)-fm(k))/(2.*eps);
      }
    }
    ipf->computeIpv(iele,0,IPState::current); // 0 for - elem and npts for + elem
    lterm->setMinus(false);
    for(int j=0;j<iele->getElem(1)->getNumVertices();j++){
      for(int jj=0;jj<3;jj++){
        if(!fullDg)
          D = Dof(iele->getElem(1)->getVertex(j)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000));
        else
          D = Dof(iele->getElem(1)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000,j));
        lterm->setPertDof(D);
        fm.scale(0.);fp.scale(0.);
        per = eps;
        ufield->set(D,per);
        lterm->setPert(per);
        ipf->computeIpv(iele,npts,IPState::current); // 0 for - elem and npts for + elem
        lterm->get(iele,npts,GP,fp);
        per = -2.*eps;
        lterm->setPert(-eps);
        ufield->set(D,per);
        ipf->computeIpv(iele,npts,IPState::current); // 0 for - elem and npts for + elem
        lterm->get(iele,npts,GP,fm);
        per = eps;
        ufield->set(D,per);
        for(int k=0;k<nbdof_m+nbdof_p;k++)
          m(k,nbdof_m+jj*iele->getElem(1)->getNumVertices()+j)+=(fp(k)-fm(k))/(2.*eps);
      }
    }
    ipf->computeIpv(iele,npts,IPState::current); // 0 for - elem and npts for + elem
    //ipf->compute1state(IPState::current);
//    m.print(" interface pert");*/
  }
  else
   printf("not implemented\n");
}

void IsotropicElasticForceInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP,fullVector<double> &m)
{
  // data initialization
  // Retrieve of the element link to the interface element velem[0] == minus elem ; velem == plus elem
  MInterfaceElement *iele = dynamic_cast<MInterfaceElement*>(ele);
  MElement ** velem = iele->getElem();
  const int nbdof_m = DgC0LinearTerm<SVector3>::space1.getNumKeys(velem[0]);
  const int nbdof_p = DgC0LinearTerm<SVector3>::space1.getNumKeys(velem[1]);
  const int nbFF_m=nbdof_m/3;
  const int nbFF_p=nbdof_p/3;
  // Initialization
  m.resize(nbdof_m+nbdof_p);
  m.scale(0.);
  double uem,uep,vem,vep;
  std::vector<TensorialTraits<double>::ValType> Val_m, Val_p;
  std::vector<TensorialTraits<double>::GradType> Grads_m;
  std::vector<TensorialTraits<double>::GradType> Grads_p;
  std::vector<TensorialTraits<double>::GradType> Grads;
  std::vector<TensorialTraits<double>::HessType> Hess_m;
  std::vector<TensorialTraits<double>::HessType> Hess_p;
  const LocalBasis* lb[3];
  double Bhat_p[256][3][2][2],Bhat_m[256][3][2][2], Bm_p[256][3][2][2],Bm_m[256][3][2][2];
  double Bs_p[256][3],Bs_m[256][3];
  LinearElasticShellHookeTensor HOOKEhat_p; LinearElasticShellHookeTensor *Hhat_p=&HOOKEhat_p;
  LinearElasticShellHookeTensor HOOKehat_m; LinearElasticShellHookeTensor *Hhat_m=&HOOKehat_m;
  LinearElasticShellHookeTensor HOOKehatmean; LinearElasticShellHookeTensor *Hhatmean=&HOOKehatmean;
  double Deltat_m[256][3][3], Deltat_p[256][3][3];
  double Cmt,Cnt,Cst,wJ;
  double me_cons[3];
  double me_comp[3];
  double me_stab[3];
  reductionElement nhatmean, mhatmean;
  SVector3 ujump;
  std::vector<bool> vbroken;
  std::vector<bool> vDeltanNeg;

  // Characteristic size of interface element
  double h_s = iele->getCharacteristicSize();

  const double Bhs = beta1/h_s;
  const double B2hs= beta2/h_s;
  const double B3hs= beta3/h_s;
  // displacement
  std::vector<double> disp;
  if(MatrixByPerturbation)
    ufield->getForPerturbation(iele,minus,pertDof,pert,disp);
  else
    ufield->get(iele,disp);

  // sum on Gauss' points
  ipf->getBroken(iele,_elemtype,vbroken,vDeltanNeg);
  //for(int i=0;i<npts;i++){ if(vbroken[i]) printf("%d true\n",iele->getNum()); else printf("%d false\n",iele->getNum());}
  for (int i = 0; i < npts; i++)
  {
    // Coordonate of Gauss' point i
    const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
    // Weight of Gauss' point i
    const double weight = GP[i].weight;
    //Compute the position (u,v) in the element of the Gauss point (to know where evaluate the shape functions)
    iele->getuvOnElem(u,uem,vem,uep,vep);
    // Compute of gradient and hessian of shape functions on interface element and on elements
    // ATTENTION after multiplication by multipliers (functionspace 276) components are in the "second line of tensor change this ??
    DgC0LinearTerm<SVector3>::space1.gradfuvw(iele,u, v, w, Grads); // grad of shape fonction on interface element
    DgC0LinearTerm<SVector3>::space1.gradfuvw(velem[0],uem, vem, w, Grads_m); // w = 0
    DgC0LinearTerm<SVector3>::space1.gradfuvw(velem[1],uep, vep, w, Grads_p);
    DgC0LinearTerm<SVector3>::space1.hessfuvw(velem[0],uem, vem, w, Hess_m);
    DgC0LinearTerm<SVector3>::space1.hessfuvw(velem[1],uep, vep, w, Hess_p);

    // m^a terms (Consistency) and get data
    if(!fullDg) ipf->getMomentReductionAndLocalBasis(iele,i,npts,_elemtype,IPState::current,mhatmean,lb);
    else ipf->getReductionAndLocalBasis(iele,i,npts,_elemtype,IPState::current,nhatmean,mhatmean,lb);

    wJ = weight*lb[2]->getJacobian();
    // Compute of Deltat_tilde
    compute_Deltat_tilde(lb[1],Grads_p,nbFF_p,Deltat_p);
    compute_Deltat_tilde(lb[0],Grads_m,nbFF_m,Deltat_m);

    // Assemblage
    double dt[3];
    for(int alpha=0;alpha<2;alpha++){
//      SVector3 malpha = mhatmean(alpha,0)*lb[2]->getphi0(0)+mhatmean(alpha,1)*lb[2]->getphi0(1);
      SVector3 malpha = mhatmean(0,alpha)*lb[2]->getphi0(0)+mhatmean(1,alpha)*lb[2]->getphi0(1);
      for(int j=0;j<nbFF_m;j++){
        matTvectprod(Deltat_m[j],malpha,dt);
        for(int k=0;k<3;k++)
          m(j+k*nbFF_m)+= - (wJ*dt[k]*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha))));
      }
      for(int j=0;j<nbFF_p;j++){
        matTvectprod(Deltat_p[j],malpha,dt);
        for(int k=0;k<3;k++)
          m(j+k*nbFF_p+nbdof_m)+= (wJ*dt[k]*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha))));
      }
    }

    // Compatibility and stability (if not broken)
    if(!vbroken[i]){
      Compute_Bm(lb[1],Grads_p,Hess_p,nbFF_p,Bm_p);
      Compute_Bm(lb[0],Grads_m,Hess_m,nbFF_m,Bm_m);
      compute_Bhat(lb[1],nbFF_p,Bm_p,Bhat_p);
      compute_Bhat(lb[0],nbFF_m,Bm_m,Bhat_m);

      // Compute of Hooke hat tensor on minus and plus element
      Cmt = wJ*Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
      Hhat_p->hat(lb[1],Cmt,nu);
      Hhat_m->hat(lb[0],Cmt,nu);
      Hhatmean->mean(Hhat_p,Hhat_m); // mean value of tensor by component used for stability term
      // Assemblage (compatibility and stability)
      for(int j=0;j<nbFF_m;j++){
        compC0PlateForceTerms(Hhat_m,nbFF_m,nbFF_p,Bhat_m[j],Deltat_m,Deltat_p,lb[2],disp,me_comp);
        stabilityC0PlateForceTerms(nbFF_m,nbFF_p,Hhatmean,Deltat_m[j],Deltat_m,Deltat_p,lb[2],disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+jj*nbFF_m) += ( me_comp[jj]- Bhs * me_stab[jj] );
      }
      for(int j=0;j<nbFF_p;j++){
        compC0PlateForceTerms(Hhat_p,nbFF_m,nbFF_p,Bhat_p[j],Deltat_m,Deltat_p,lb[2],disp,me_comp);
        stabilityC0PlateForceTerms(nbFF_m,nbFF_p,Hhatmean,Deltat_p[j],Deltat_m,Deltat_p,lb[2],disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+(jj*nbFF_p)+nbdof_m) += ( me_comp[jj] + Bhs * me_stab[jj]);
      }
    }
    // Add extra terms for fullDg formulation
    if(fullDg){
      // Shape functions evaluated in u,v
      DgC0LinearTerm<SVector3>::space1.fuvw(velem[0],uem, vem, w, Val_m);
      DgC0LinearTerm<SVector3>::space1.fuvw(velem[1],uep, vep, w, Val_p);

      // n^a terms
      //ipf->getStressReduction(iele,i,npts,_elemtype,IPState::current,nhatmean); did before
      // Assemblage (Consistency)
      for(int alpha=0;alpha<2;alpha++){
        SVector3 nalpha = nhatmean(alpha,0) * lb[2]->getphi0(0) + nhatmean(alpha,1)*lb[2]->getphi0(1);
        for(int j=0;j<nbFF_m;j++)
          for(int k=0;k<3;k++)
            m(j+k*nbFF_m)+=- (wJ*nalpha(k)*Val_m[j]*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha))));
        for(int j=0;j<nbFF_p;j++)
          for(int k=0;k<3;k++)
            m(j+k*nbFF_p+nbdof_m)+=(wJ*nalpha(k)*Val_p[j]*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha))));
      }

      // Compatibility and stability (if not broken)
      if(!vbroken[i]){
        // compute jump of u
        displacementjump(Val_m,nbFF_m,Val_p,nbFF_p,disp,ujump);
        // Hooke tensor (symmetrization is for elastic part)
        Cnt = wJ*Cn;
        Hhat_m->set(lb[0],Cnt,nu);
        Hhat_p->set(lb[1],Cnt,nu);
        Hhat_m->hat(lb[0],Cnt,nu);
        Hhat_p->hat(lb[1],Cnt,nu);
        Hhatmean->mean(Hhat_m,Hhat_p);
        // Assemblage (Compatibility and stability)
        for(int alpha=0;alpha<2;alpha++){
          for(int beta=0;beta<2;beta++)
            for(int gamma=0;gamma<2;gamma++)
              for(int delta=0;delta<2;delta++){
                stabC0PlateForceMembraneTerms(beta,delta,lb[2],ujump,me_stab);
                double tempStab = Hhatmean->get(alpha,beta,gamma,delta)*B2hs*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha)))*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(gamma)));
                double tempCompm = 0.5*Hhat_m->get(alpha,beta,gamma,delta)*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha)));
                double tempCompp = 0.5*Hhat_p->get(alpha,beta,gamma,delta)*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha)));
                for(int j=0;j<nbFF_m;j++){
                  compC0PlateForceMembraneTerms(beta,gamma,delta,lb[0],lb[2],Grads_m[j],ujump,me_comp);
                  for(int jj=0;jj<3;jj++)
                    m(j+jj*nbFF_m)+=(tempCompm*me_comp[jj]-tempStab*Val_m[j]*me_stab[jj]);
                }
                for(int j=0;j<nbFF_p;j++){
                  compC0PlateForceMembraneTerms(beta,gamma,delta,lb[1],lb[2],Grads_p[j],ujump,me_comp);
                  for(int jj=0;jj<3;jj++)
                    m(j+jj*nbFF_p+nbdof_m)+=(tempCompp*me_comp[jj]+tempStab*Val_p[j]*me_stab[jj]);
                }
              }
        }
      // Out of Plane term (if not totally broken)
      Cst = wJ*Cs*B3hs; // Hooke tensor in "shearing" 1 component = Cs
      // compute B vector
      compute_Bs(lb[2],Val_m,nbFF_m,Bs_m);
      compute_Bs(lb[2],Val_p,nbFF_p,Bs_p);
      for(int j=0;j<nbFF_m;j++){
        stabilityC0PlateForceShearingTerms(Bs_m[j], Bs_m, Bs_p,nbFF_m,nbFF_p,disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+jj*nbFF_m) += -Cst * me_stab[jj];
      }
      for(int j=0;j<nbFF_p;j++){
        stabilityC0PlateForceShearingTerms(Bs_p[j],Bs_m,Bs_p,nbFF_m,nbFF_p,disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+(jj*nbFF_p)+nbdof_m) += Cst * me_stab[jj];
      }
    }
    Val_m.clear(); Val_p.clear();
  }
  // Because component are push_back in Grads in gradfuvw idem for hess
  Grads_m.clear(); Grads_p.clear(); Hess_m.clear(); Hess_p.clear(); Grads.clear();
  }
  if(inverseSign)
    for(int i=0;i<nbdof_m+nbdof_p;i++) m(i) = -m(i);
}

void IsotropicElasticStiffVirtualInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m)
 {
  if (sym)
  {
    // data initialization
    // Retrieve of the element link to the interface element velem[0] == minus elem ; velem == plus elem
    MInterfaceElement *iele = dynamic_cast<MInterfaceElement*>(ele);
    MElement ** velem = iele->getElem();
    const int nbdof_m = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(velem[0]);
    const int nbFF_m = nbdof_m/3;
    // Initialization
    m.resize(nbdof_m, nbdof_m);
    m.setAll(0.);
    double uem,uep,vem,vep;
    std::vector<TensorialTraits<double>::GradType> Grads_m;
    std::vector<TensorialTraits<double>::GradType> Grads_p;
    std::vector<TensorialTraits<double>::GradType> Grads;
    std::vector<TensorialTraits<double>::HessType> Hess_m;
    std::vector<TensorialTraits<double>::HessType> Hess_p;
    LocalBasis LBS,LBM;
    LocalBasis *lbs=&LBS; LocalBasis *lb_m=&LBM;
    //std::vector<std::vector<fullMatrix<double> > > Bhat_m;
    double Bhat_m[256][3][2][2], Bm_m[256][3][2][2];
    LinearElasticShellHookeTensor HOOKehat_m;
    LinearElasticShellHookeTensor *Hhat_m=&HOOKehat_m;
    double Deltat_m[256][3][3],Deltat_p[256][3][3];
    double Cmt;
    double me_cons[3][3];
    double me_comp[3][3];
    double me_stab[3][3];

    // Characteristic size of interface element
    double h_s = iele->getCharacteristicSize();
    const double Bhs = beta1/h_s;
    // sum on Gauss' points
    for (int i = 0; i < npts; i++)
    {
      // Coordonate of Gauss' point i
      const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
      // Weight of Gauss' point i
      const double weight = GP[i].weight;

      //printf("Abscisse of gauss point %f\n",u);
      //Compute the position (u,v) in the element of the Gauss point (to know where evaluate the shape functions)
      iele->getuvOnElem(u,uem,vem,uep,vep);
      //printf("Position (u,v) of minus element (%f,%f)\n",uem,vem);
      // Compute of gradient and hessian of shape functions on interface element and on elements
      // ATTENTION after multiplication by multipliers (functionspace 276) components are in the "second line of tensor change this ??
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(iele,u, v, w, Grads); // grad of shape fonction on interface element
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(velem[0],uem, vem, w, Grads_m); // w = 0
      DgC0BilinearTerm<SVector3,SVector3>::space1.hessfuvw(velem[0],uem, vem, w, Hess_m);

      // basis of elements and interface element
      lb_m->set(velem[0],Grads_m,Hess_m); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lbs->set(iele,Grads,lb_m->gett0()); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??

      // PushForwardTensor
      lb_m->set_pushForward(lbs);

      // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
      Compute_Bm(lb_m,Grads_m,Hess_m,nbFF_m,Bm_m);
      compute_Bhat(lb_m,nbFF_m,Bm_m,Bhat_m);
      // Compute of Hooke hat tensor on minus and plus element
      Cmt = weight * lbs->getJacobian()  * Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
      Hhat_m->hat(lb_m,Cmt,nu);

      // Compute of Deltat_tilde
      compute_Deltat_tildeBound(lb_m,Grads_m,nbFF_m,Deltat_m,lbs);
      // loop on dof ATTENTION SAME NUMBER OF DOF for the two elements TODO take into account a different dof numbers between the two elements There ok because sym ??
      for(int j=0;j<3;j++) for(int jj=0;jj<3;jj++) me_comp[j][jj] = me_stab[j][jj]=0.;
      for(int j=0;j<nbFF_m;j++){
        for(int k=0;k<nbFF_m;k++){
          consC0PlateStiffnessTerms(Hhat_m,Bhat_m[k],Deltat_m[j],lbs,me_cons); // Error Fix IT
          compC0PlateStiffnessTerms(Hhat_m,Bhat_m[j],Deltat_m[k],lbs,me_comp);
          stabilityC0PlateStiffnessTerms(Hhat_m,Deltat_m[j], Deltat_m[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+jj*nbFF_m,k+kk*nbFF_m) += -(me_cons[jj][kk] + me_comp[jj][kk] - Bhs * me_stab[jj][kk]);
        }
      }

      // Because component are push_back in Grads in gradfuvw idem for hess
      Grads_m.clear(); Hess_m.clear(); Grads.clear();
  }

/*   m.print("Virtual InterfaceBound");
   m.setAll(0.);
    // By numerical perturbation (Verification OK)
    double eps=1.e-5;
    double per;
    fullVector<double> fp(nbdof_m);
    fullVector<double> fm(nbdof_m);
    fp.scale(0.);
    fm.scale(0.);
    DgC0LinearTerm<SVector3> *lterm = this->getLinearTerm();
    Dof D(0,0);
    for(int j=0;j<iele->getElem(0)->getNumVertices();j++){
      for(int jj=0;jj<3;jj++){
        if(!fullDg)
         D = Dof(iele->getElem(0)->getVertex(j)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000));
        else
         D = Dof(iele->getElem(0)->getNum(),DgC0PlateDof::createTypeWithThreeInts(jj,1000,j));
        fm.scale(0.);fp.scale(0.);
        per = eps;
        ufield->set(D,per);
        ipf->compute1state(IPState::current);
        lterm->get(iele,npts,GP,fp);
        per = -2.*eps;
        ufield->set(D,per);
        ipf->compute1state(IPState::current);
        lterm->get(iele,npts,GP,fm);
        per = eps;
        ufield->set(D,per);
        for(int k=0;k<nbdof_m;k++)
          m(k,jj*iele->getElem(0)->getNumVertices()+j)=(fp(k)-fm(k))/(2.*eps);
      }
    }
    ipf->compute1state(IPState::current);
    m.print("Virtual interface pert");*/
  }
  else
    printf("not implemented\n");
}

void IsotropicElasticForceVirtualInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP, fullVector<double> &m)
{
  // data initialization
  // Retrieve of the element link to the interface element velem[0] == minus elem ; velem == plus elem
  MInterfaceElement *iele = dynamic_cast<MInterfaceElement*>(ele);
  MElement ** velem = iele->getElem();
  const int nbdof_m = DgC0LinearTerm<SVector3>::space1.getNumKeys(velem[0]);
  const int nbFF_m = nbdof_m/3;
  // Initialization
  m.resize(nbdof_m);
  m.scale(0.);
  double uem,vem,uep,vep;
  std::vector<TensorialTraits<double>::GradType> Grads_m;
  std::vector<TensorialTraits<double>::GradType> Grads;
  std::vector<TensorialTraits<double>::HessType> Hess_m;
  const LocalBasis *lb[3];
  double Bhat_m[256][3][2][2], Bm_m[256][3][2][2];
  LinearElasticShellHookeTensor HOOKehat_m;
  LinearElasticShellHookeTensor *Hhat_m=&HOOKehat_m;
  double Deltat_m[256][3][3];
  double Cmt,wJ;
  double me_cons[3];
  double me_comp[3];
  double me_stab[3];
  reductionElement mhatmean;
  // Characteristic size of interface element
  double h_s = iele->getCharacteristicSize();

  const double Bhs = beta1/h_s;

  // displacement
  std::vector<double> disp;
  ufield->get(iele,disp);
  // sum on Gauss' points
  for(int i = 0; i < npts; i++)
  {
    // Coordonate of Gauss' point i
   const double u = GP[i].pt[0]; const double v = GP[i].pt[1]; const double w = GP[i].pt[2];
    // Weight of Gauss' point i
    const double weight = GP[i].weight;
    //printf("Abscisse of gauss point %f\n",u);
    //Compute the position (u,v) in the element of the Gauss point (to know where evaluate the shape functions)
    iele->getuvOnElem(u,uem,vem,uep,vep);
    // Compute of gradient and hessian of shape functions on interface element and on elements
    // ATTENTION after multiplication by multipliers (functionspace 276) components are in the "second line of tensor change this ??
    DgC0LinearTerm<SVector3>::space1.gradfuvw(iele,u, v, w, Grads); // grad of shape fonction on interface element
    DgC0LinearTerm<SVector3>::space1.gradfuvw(velem[0],uem, vem, w, Grads_m); // w = 0
    DgC0LinearTerm<SVector3>::space1.hessfuvw(velem[0],uem, vem, w, Hess_m);


     // get data
    ipf->getVirtualMomentReductionAndLocalBasis(iele,i,npts,_elemtype,IPState::current,mhatmean,lb);
     wJ=weight * lb[2]->getJacobian();
    // Compute of Deltat_tilde
    compute_Deltat_tildeBound(lb[0],Grads_m,nbFF_m,Deltat_m,lb[2]);

    // Consistency
    double dtmalpha[3];
    for(int alpha=0;alpha<2;alpha++){
      SVector3 malpha = mhatmean(alpha,0)*lb[2]->getphi0(0)+mhatmean(alpha,1)*lb[2]->getphi0(1);
      for(int j=0;j<nbFF_m;j++){
        matTvectprod(Deltat_m[j],malpha,dtmalpha);
        for(int k=0;k<3;k++)
          m(j+k*nbFF_m)+= - (wJ*dtmalpha[k]*(-scaldot(lb[2]->getphi0(1),lb[2]->getphi0(alpha))));
      }
    }

    // Compatibility and stability
    // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
    Compute_Bm(lb[0],Grads_m,Hess_m,nbFF_m,Bm_m);
    compute_Bhat(lb[0],nbFF_m,Bm_m,Bhat_m);
    Cmt = wJ * Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
    Hhat_m->hat(lb[0],Cmt,nu);
    // Assemblage
    for(int j=0;j<nbFF_m;j++){
      compC0PlateForceTermsBound(Hhat_m,nbFF_m,Bhat_m[j],Deltat_m,lb[2],disp,me_comp);
      stabilityC0PlateForceTermsBound(nbFF_m,Hhat_m,Deltat_m[j],Deltat_m,lb[2],disp,me_stab);
      for(int jj=0;jj<3;jj++){
        m(j+jj*nbFF_m) += (me_comp[jj] - Bhs * me_stab[jj] );
      }
    }
    // Because component are push_back in Grads in gradfuvw idem for hess
    Grads_m.clear(); Hess_m.clear(); Grads.clear();
  }
  if(inverseSign)
    for(int i=0;i<nbdof_m;i++) m(i) =-m(i);
};
