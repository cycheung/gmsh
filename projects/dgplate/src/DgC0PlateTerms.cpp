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
#include "DgC0PlateElementaryTerms.h"
// the following functions are defined in DgC0PlateElementaryTerms.h
inline void diaprod(const double a[3], const double b[3], double m[3][3]);
inline double scaldot(const double a[3],const SVector3 b);
static inline void dot(const double a[3], const double b[3], double c[3]);
static inline void dot(const double a[3], const SVector3 &b, double c[3]);
inline void matTvectprod(const double m[3][3], const SVector3 &v, double v2[3]);
inline void matTvectprod(const double m[3][3], const double v[3],  double v2[3]);
static inline void matvectprod(const double m[3][3], const SVector3 &v1,double v2[3]);
inline void BulkC0PlateDGStiffnessMembraneTerms(const double Bj[3][2][2],const double Bk[3][2][2], const LinearElasticShellHookeTensor *H, double me[3][3]);
inline void BulkC0PlateDGStiffnessBendingTerms(TensorialTraits<double>::HessType &hessj, TensorialTraits<double>::HessType &hessk, const LinearElasticShellHookeTensor *H, const LocalBasis *lb, double me[3][3]);
inline void consC0PlateStiffnessTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double dt[3][3], const LocalBasis *lb, double me[3][3]);
inline void compC0PlateStiffnessTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double dt[3][3], const LocalBasis *lb, double me[3][3]);
inline void stabilityC0PlateStiffnessTerms(LinearElasticShellHookeTensor *Hhat, const double dta[3][3], const double dtb[3][3], const LocalBasis *lb, double me[3][3]);
inline void consC0PlateStiffnessMembraneTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double Na, const LocalBasis *lb, double me[3][3]);
inline void compC0PlateStiffnessMembraneTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const double Na, const LocalBasis *lb, double me[3][3]);
inline void stabilityC0PlateStiffnessMembraneTerms(LinearElasticShellHookeTensor *Hhat, const double Na, const double Nb, const LocalBasis *lb, double me[3][3]);
inline void stabilityC0PlateStiffnessShearingTerms(const double Bj[3], const double Bk[3],double me[3][3]);
inline void consC0PlateForceMembraneTerms(LinearElasticShellHookeTensor *Hhat,const double Bhat[3][2][2],const std::vector<SVector3> &Vals_m,const std::vector<SVector3> &Vals_p, const LocalBasis *lb, const fullMatrix<double> &disp,double me[3]);
inline void compC0PlateForceMembraneTerms(LinearElasticShellHookeTensor *Hhat_m,const LinearElasticShellHookeTensor *Hhat_p,const double Bhat_m[256][3][2][2],const double Bhat_p[256][3][2][2], const int n_m,const int n_p, const double Na, const LocalBasis *lb,const fullMatrix<double> &disp, double me[3]);
inline void compC0PlateForceMembraneTerms(const int beta,const int gamma,const int delta,const LocalBasis *lb, const LocalBasis *lbs,const TensorialTraits<double>::GradType &Gradj,const SVector3 &ujump,double me_comp[3]);
inline void stabilityC0PlateForceMembraneTerms(LinearElasticShellHookeTensor *Hhat, const std::vector<SVector3> &Vals_m,const std::vector<SVector3> &Vals_p, const double Nb, const LocalBasis *lb, const fullMatrix<double> &disp,double me[3]);
inline void stabC0PlateForceMembraneTerms(const int beta,const int delta,const LocalBasis *lb,const SVector3 &ujump,double me[3]);
inline void stabilityC0PlateForceShearingTerms(const double Bj[3], const double B_m[256][3],const double B_p[256][3],const int n_m,const int n_p,const fullMatrix<double> &disp,double me[3]);
inline void consC0PlateForceTerms(const LinearElasticShellHookeTensor *Hhat, const int n_m, const int n_p, const double Bhat[3][2][2], const double Dt_m[256][3][3], const double Dt_p[256][3][3],const LocalBasis *lb, const fullMatrix<double> &disp, double me[3]);
inline void compC0PlateForceTerms(const int n_m, const int n_p, const LinearElasticShellHookeTensor *Hhat_m, const LinearElasticShellHookeTensor *Hhat_p, const double Bhat_m[256][3][2][2],const double Bhat_p[256][3][2][2],const double Dt[3][3],const LocalBasis *lb, const fullMatrix<double> &disp, double me[3]);
inline void stabilityC0PlateForceTerms(const int n_p,const int n_m, const LinearElasticShellHookeTensor *Hhat,const  double Dt[3][3],const double Dt_m[256][3][3],const double Dt_p[256][3][3], const LocalBasis *lb, const fullMatrix<double> &disp, double me[3]);
inline void consC0PlateForceTermsBound(const LinearElasticShellHookeTensor *Hhat, const int n, const double Bhat[3][2][2], const double Dt_m[256][3][3], const LocalBasis *lb, const fullMatrix<double> &disp, double me[3]);
inline void compC0PlateForceTermsBound(const int n_m,const LinearElasticShellHookeTensor *Hhat_m, const double Bhat_m[256][3][2][2], const double Dt[3][3],const LocalBasis *lb, const fullMatrix<double> &disp, double me[3]);
inline void stabilityC0PlateForceTermsBound(const int n_m,const LinearElasticShellHookeTensor *Hhat,const double Dt[3][3],const double Dt_m[256][3][3], const LocalBasis *lb, const fullMatrix<double> &disp, double me[3]);
void compute_Bn(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int n, double B[][3][2][2]);
void compute_Bs(const LocalBasis *lb, const std::vector<TensorialTraits<double>::ValType> &Vals, const int n, double B[][3]);
void compute_Bnhat(const LocalBasis *lb, const int n, const double B[256][3][2][2], double Bhat[256][3][2][2]);
void  Compute_Bhat(const LocalBasis *lb,const std::vector<TensorialTraits<double>::HessType> &Hess, const int &n, double B[256][3][2][2]);
void compute_Deltat_tilde(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int &n, double Deltat[256][3][3]);
void compute_Deltat_tildeBound(const LocalBasis *lb, const std::vector<TensorialTraits<double>::GradType> &Grads, const int &n, double Deltat[256][3][3], const LocalBasis *lbs);


void IsotropicElasticBulkTermC0Plate::get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m)
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
    double Bn[256][3][2][2]; // max order 256 or dynamical allocation ?? better than std::vector and fullMatrix ?? create a class with this variable
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

      lb->set(ele,Grad); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> qui retourne un vecteur vide prob de template??

      // multiplication of constant by detJ and weight
      Cmt = lb->getJacobian() * weight *Cm;
      Cnt = lb->getJacobian() * weight *Cn;
      // compute of Hooke tensor
      H->set(lb,nu);

      // compute Bn value
      compute_Bn(lb,Grad,nbFF,Bn);

     // loop on SF to construct the elementary (bulk) stiffness matrix at the Gauss' point i
     for(int j=0; j<nbFF;j++)
       for(int k=0;k<nbFF;k++){
         BulkC0PlateDGStiffnessBendingTerms(Hess[j],Hess[k],H,lb,me_bending);
         BulkC0PlateDGStiffnessMembraneTerms(Bn[j],Bn[k],H,me_membrane);
         for(int jj=0;jj<3;jj++)
           for(int kk=0;kk<3;kk++)
             m(j+(jj*nbFF),k+(kk*nbFF)) += (Cmt*me_bending[jj][kk]+Cnt*me_membrane[jj][kk]);
       }
     // clear the hessian and Grad because the components append in hessfuvw and gradfuvw
     Hess.clear(); Grad.clear();
    }
/*    m.print("bulk");
    // By numerical perturbation (Verification OK)
    double eps=1.e-8;
    fullMatrix<double> fp(nbdof,1);
    fullMatrix<double> fm(nbdof,1);
    fp.setAll(0.);
    fm.setAll(0.);
    fullMatrix<double> m2(nbdof,nbdof); m2.setAll(0.);
    fullMatrix<double> disp(nbdof,1);
    for(int j=0;j<nbdof;j++){
      disp.setAll(0.); fm.setAll(0.);fp.setAll(0.);
      disp(j,0)=eps;
      this->getForce(ele,npts,GP,disp,fp);
      disp(j,0)=-eps;
      this->getForce(ele,npts,GP,disp,fm);
      for(int k=0;k<nbdof;k++)
        m2(k,j)=(fp(k,0)-fm(k,0))/(2.*eps);
    }
    m2.print("Bulk pert");*/
  }
  else
    printf("not implemented\n");
}

void IsotropicElasticBulkTermC0Plate::get(MElement *ele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m)
{
  if (sym)
  {
    // Initialization of some data
    const int nbdof = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(ele);
    const int nbFF = nbdof/3;
    double Cmt= 0.,Cnt=0.;
    LinearElasticShellHookeTensor HOOKe; LinearElasticShellHookeTensor *H=&HOOKe; // make a pointer in an other way
    m.resize(nbdof,1);
    m.setAll(0.);
    std::vector<TensorialTraits<double>::HessType> Hess;
    std::vector<TensorialTraits<double>::GradType> Grad;
    LocalBasis LB;
    LocalBasis *lb=&LB; // two last line in 1 operation ??
    std::vector<SVector3> nalpha,malpha;
    nalpha.reserve(2); malpha.reserve(2);

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

      lb->set(ele,Grad); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> qui retourne un vecteur vide prob de template??

      // multiplication of constant by detJ and weight
      Cmt = lb->getJacobian() * weight * Cm;
      Cnt = lb->getJacobian() * weight * Cn;

      // compute of Hooke tensor
      H->set(lb,nu);

      stressReduction(H,Grad,lb,disp,0,nalpha);
      momentReduction(H,Hess,lb,disp,0,malpha);
      for(int j=0; j<nbFF;j++){
        for(int k=0;k<3;k++)
          for(int alpha=0;alpha<2;alpha++)
            for(int beta=0;beta<2;beta++)
              m(j+k*nbFF,0) += (Cnt*Grad[j](alpha)*nalpha[alpha](beta)*lb->getphi0(beta,k)- Cmt*Hess[j](alpha,beta)*malpha[alpha](beta)*lb->gett0(k));
      }
    // clear the hessian and Grad because the components append in hessfuvw and gradfuvw
    Hess.clear(); Grad.clear();
    }
  }
  else
    printf("not implemented\n");
}

void IsotropicElasticInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP, fullMatrix<double> &m)
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
    double Bhat_p[256][3][2][2], Bhat_m[256][3][2][2], Bn_m[256][3][2][2], Bn_p[256][3][2][2];
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
      lb_m->set(velem[0],Grads_m); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lb_p->set(velem[1],Grads_p); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lbs->set(iele,Grads,lb_p->gett0(),lb_m->gett0()); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      // PushForwardTensor
      lb_m->set_pushForward(lbs);
      lb_p->set_pushForward(lbs);

      // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
      Compute_Bhat(lb_p,Hess_p,nbFF_p,Bhat_p);
      Compute_Bhat(lb_m,Hess_m,nbFF_m,Bhat_m);

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
          consC0PlateStiffnessTerms(Hhat_m,Bhat_m[j],Deltat_m[k],lbs,me_cons);
          compC0PlateStiffnessTerms(Hhat_m,Bhat_m[k],Deltat_m[j],lbs,me_comp);
          stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_m[j], Deltat_m[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+jj*nbFF_m,k+kk*nbFF_m) += (- me_cons[jj][kk] - me_comp[jj][kk] + Bhs * me_stab[jj][kk] );
        }
        for(int k=0;k<nbFF_p;k++){
          consC0PlateStiffnessTerms(Hhat_m,Bhat_m[j],Deltat_p[k],lbs,me_cons);
          compC0PlateStiffnessTerms(Hhat_p,Bhat_p[k],Deltat_m[j],lbs,me_comp);
          stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_m[j], Deltat_p[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+jj*nbFF_m,k+nbdof_m+kk*nbFF_p) += ( me_cons[jj][kk] - me_comp[jj][kk]- Bhs * me_stab[jj][kk] );
        }
      }
      for(int j=0;j<nbFF_p;j++){
        for(int k=0;k<nbFF_m;k++){
          consC0PlateStiffnessTerms(Hhat_p,Bhat_p[j],Deltat_m[k],lbs,me_cons);
          compC0PlateStiffnessTerms(Hhat_m,Bhat_m[k],Deltat_p[j],lbs,me_comp);
          stabilityC0PlateStiffnessTerms(Hhatmean,Deltat_p[j], Deltat_m[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+(jj*nbFF_p)+nbdof_m,k+kk*nbFF_m) += (- me_cons[jj][kk] + me_comp[jj][kk] - Bhs * me_stab[jj][kk]);
        }
        for(int k=0;k<nbFF_p;k++){
          consC0PlateStiffnessTerms(Hhat_p,Bhat_p[j],Deltat_p[k],lbs,me_cons);
          compC0PlateStiffnessTerms(Hhat_p,Bhat_p[k],Deltat_p[j],lbs,me_comp);
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
      compute_Bnhat(lb_m,nbFF_m,Bn_m,Bhat_m);
      compute_Bnhat(lb_p,nbFF_p,Bn_p,Bhat_p);

      // Compute Hooke tensor
      Cnt = weight*lbs->getJacobian()*Cn;
      Hhat_m->hat(lb_m,Cnt,nu); // Redondant avec + haut (faut enlever le Cmt pour le calcul augmente cout calcul ??)
      Hhat_p->hat(lb_p,Cnt,nu);
      Hhatmean->mean(Hhat_m,Hhat_p);
      for(int j=0;j<3;j++) for(int k=0;k<3;k++){me_stab[j][k]=0.;me_cons[j][k]=0.;}

      // Assembly (regroup all ??)
      for(int j=0;j<nbFF_m;j++){
        for(int k=0;k<nbFF_m;k++){
          consC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_m[j],Val_m[k],lbs,me_cons);
          compC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_m[k],Val_m[j],lbs,me_comp);
          stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_m[j], Val_m[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+jj*nbFF_m,k+kk*nbFF_m) += (- me_cons[jj][kk] - me_comp[jj][kk] + B2hs * me_stab[jj][kk] );
        }
        for(int k=0;k<nbFF_p;k++){
          consC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_m[j],Val_p[k],lbs,me_cons);
          compC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_p[k],Val_m[j],lbs,me_comp);
          stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_m[j], Val_p[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+jj*nbFF_m,k+nbdof_m+kk*nbFF_p) += ( me_cons[jj][kk] - me_comp[jj][kk]- B2hs * me_stab[jj][kk] );
        }
      }
      for(int j=0;j<nbFF_p;j++){
        for(int k=0;k<nbFF_m;k++){
          consC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_p[j],Val_m[k],lbs,me_cons);
          compC0PlateStiffnessMembraneTerms(Hhat_m,Bhat_m[k],Val_p[j],lbs,me_comp);
          stabilityC0PlateStiffnessMembraneTerms(Hhatmean,Val_p[j], Val_m[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+(jj*nbFF_p)+nbdof_m,k+kk*nbFF_m) += (- me_cons[jj][kk] + me_comp[jj][kk] - B2hs * me_stab[jj][kk]);
        }
        for(int k=0;k<nbFF_p;k++){
          consC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_p[j],Val_p[k],lbs,me_cons);
          compC0PlateStiffnessMembraneTerms(Hhat_p,Bhat_p[k],Val_p[j],lbs,me_comp);
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
              {m(j+jj*nbFF_m,k+kk*nbFF_m) += Cst * me_stab[jj][kk]; }//printf("-- %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);}
        }
        for(int k=0;k<nbFF_p;k++){
          stabilityC0PlateStiffnessShearingTerms(Bs_m[j], Bs_p[k],me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              {m(j+jj*nbFF_m,k+nbdof_m+kk*nbFF_p) += -Cst * me_stab[jj][kk]; }//printf("-+ %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);}
        }
      }
      for(int j=0;j<nbFF_p;j++){
        for(int k=0;k<nbFF_m;k++){
          stabilityC0PlateStiffnessShearingTerms(Bs_p[j], Bs_m[k],me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              {m(j+(jj*nbFF_p)+nbdof_m,k+kk*nbFF_m) += -Cst * me_stab[jj][kk]; }//printf("+- %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);}
        }
        for(int k=0;k<nbFF_p;k++){
          stabilityC0PlateStiffnessShearingTerms(Bs_p[j], Bs_p[k],me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              {m(j+(jj*nbFF_p)+nbdof_m,k+(kk*nbFF_p)+nbdof_m) += Cst * me_stab[jj][kk]; }//printf("++ %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);}
        }
      }
      Val_m.clear(); Val_p.clear();
    }

    // Because component are push_back in Grads in gradfuvw idem for hess
    Grads_m.clear(); Grads_p.clear(); Hess_m.clear(); Hess_p.clear(); Grads.clear();
  }
/*  m.print("Interface");
  // Compute the matrix by numerical perturbation Verif OK
  double eps=1.e-8;
  fullMatrix<double> fp(nbdof_m+nbdof_p,1);
  fullMatrix<double> fm(nbdof_m+nbdof_p,1);
  fp.setAll(0.);
  fm.setAll(0.);
  fullMatrix<double> m2(nbdof_m+nbdof_p,nbdof_m+nbdof_p); m2.setAll(0.);
  fullMatrix<double> disp(nbdof_m+nbdof_p,1);
  for(int j=0;j<nbdof_m+nbdof_p;j++){
    disp.setAll(0.); fm.setAll(0.);fp.setAll(0.);
    disp(j,0)=eps;
    this->getInterForce(iele,npts,GP,disp,fullDg,fp);
    disp(j,0)=-eps;
    this->getInterForce(iele,npts,GP,disp,fullDg,fm);
    for(int k=0;k<nbdof_m+nbdof_p;k++)
      m2(k,j)=(fp(k,0)-fm(k,0))/(2.*eps);
  }
  m2.print("Matrix by perturbation");*/
  /*fullMatrix<double> m3(nbdof_m+nbdof_p,nbdof_m+nbdof_p);
  for(int ii=0;ii<nbdof_m+nbdof_p;ii++)
    for(int jj=0;jj<nbdof_m+nbdof_p;jj++)
  m3(ii,jj)=m2(ii,jj)-m(ii,jj);
  m3.print("diff");*/
  /*m=m2;*/
  }
  else
   printf("not implemented\n");
}

void IsotropicElasticInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP,const fullMatrix<double> &disp,fullMatrix<double> &m)
{
  if (sym)
  {
    // data initialization
    // Retrieve of the element link to the interface element velem[0] == minus elem ; velem == plus elem
    MInterfaceElement *iele = dynamic_cast<MInterfaceElement*>(ele);
    MElement ** velem = iele->getElem();
    const int nbdof_m = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(velem[0]);
    const int nbdof_p = DgC0BilinearTerm<SVector3,SVector3>::space1.getNumKeys(velem[1]);
    const int nbFF_m=nbdof_m/3;
    const int nbFF_p=nbdof_p/3;
    // Initialization
    m.resize(nbdof_m+nbdof_p, 1);
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
    double Bhat_p[256][3][2][2],Bhat_m[256][3][2][2],Bn_m[256][3][2][2], Bn_p[256][3][2][2];
    double Bs_p[256][3],Bs_m[256][3];
    LinearElasticShellHookeTensor HOOKEhat_p; LinearElasticShellHookeTensor *Hhat_p=&HOOKEhat_p;
    LinearElasticShellHookeTensor HOOKehat_m; LinearElasticShellHookeTensor *Hhat_m=&HOOKehat_m;
    LinearElasticShellHookeTensor HOOKehatmean; LinearElasticShellHookeTensor *Hhatmean=&HOOKehatmean;
    double Deltat_m[256][3][3], Deltat_p[256][3][3];
    double Cmt,Cnt,Cst;
    double me_cons[3];
    double me_comp[3];
    double me_stab[3];
    std::vector<SVector3> nalpha_m,nalpha_p;
    nalpha_m.reserve(2);
    //nalpha_m.push_back(SVector3(0.,0.,0.));
    //nalpha_m.push_back(SVector3(0.,0.,0.));
    nalpha_p.reserve(2);
    //nalpha_p.push_back(SVector3(0.,0.,0.));
    //nalpha_p.push_back(SVector3(0.,0.,0.));
    std::vector<SVector3> nhat_m,nhat_p;
    nhat_m.reserve(2);nhat_p.reserve(2);
    SVector3 ujump;

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
      lb_m->set(velem[0],Grads_m); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lb_p->set(velem[1],Grads_p); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lbs->set(iele,Grads,lb_p->gett0(),lb_m->gett0()); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??

      // PushForwardTensor
      lb_m->set_pushForward(lbs);
      lb_p->set_pushForward(lbs);

      // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
      Compute_Bhat(lb_p,Hess_p,nbFF_p,Bhat_p);
      Compute_Bhat(lb_m,Hess_m,nbFF_m,Bhat_m);
      // Compute of Hooke hat tensor on minus and plus element
      Cmt = weight * lbs->getJacobian()  * Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
      Hhat_p->hat(lb_p,Cmt,nu);
      Hhat_m->hat(lb_m,Cmt,nu);
      Hhatmean->mean(Hhat_p,Hhat_m); // mean value of tensor by component used for stability term

      // Compute of Deltat_tilde
      compute_Deltat_tilde(lb_p,Grads_p,nbFF_p,Deltat_p);
      compute_Deltat_tilde(lb_m,Grads_m,nbFF_m,Deltat_m);
      for(int i=0;i<3;i++) me_stab[i]=0.;

      for(int j=0;j<nbFF_m;j++){
        consC0PlateForceTerms(Hhat_m,nbFF_m,nbFF_p,Bhat_m[j],Deltat_m,Deltat_p,lbs,disp,me_cons);
        compC0PlateForceTerms(nbFF_m,nbFF_p, Hhat_m,Hhat_p,Bhat_m,Bhat_p,Deltat_m[j],lbs,disp,me_comp);
        stabilityC0PlateForceTerms(nbFF_m,nbFF_p,Hhatmean,Deltat_m[j],Deltat_m,Deltat_p,lbs,disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+jj*nbFF_m,0) += (me_cons[jj] - me_comp[jj]- Bhs * me_stab[jj] );
      }
      for(int j=0;j<nbFF_p;j++){
        consC0PlateForceTerms(Hhat_p,nbFF_m,nbFF_p,Bhat_p[j],Deltat_m,Deltat_p,lbs,disp,me_cons);
        compC0PlateForceTerms(nbFF_m,nbFF_p,Hhat_m,Hhat_p,Bhat_m,Bhat_p,Deltat_p[j],lbs,disp,me_comp);
        stabilityC0PlateForceTerms(nbFF_m,nbFF_p,Hhatmean,Deltat_p[j],Deltat_m,Deltat_p,lbs,disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+(jj*nbFF_p)+nbdof_m,0) += ( me_cons[jj] + me_comp[jj] + Bhs * me_stab[jj]);
      }
      // Add extra terms for fullDg formulation
      if(fullDg){
        // Shape functions evaluated in u,v
        DgC0BilinearTerm<SVector3,SVector3>::space1.fuvw(velem[0],uem, vem, w, Val_m);
        DgC0BilinearTerm<SVector3,SVector3>::space1.fuvw(velem[1],uep, vep, w, Val_p);

        // n^a terms
        Cnt = weight*lbs->getJacobian()*Cn; //TODO include j in computation of n ??
        Hhat_m->set(lb_m,Cnt,nu);
        Hhat_p->set(lb_p,Cnt,nu);

        stressReduction(Hhat_m,Grads_m,lb_m,disp,0,nalpha_m); // Ok verif avec bulk term perturbation matrix
        stressReduction(Hhat_p,Grads_p,lb_p,disp,nbdof_m,nalpha_p); // Ok verif avec bulk term perturbation matrix
        // Rotation
        stressReductionHat(nalpha_m,lb_m,nhat_m);
        stressReductionHat(nalpha_p,lb_p,nhat_p);

        // Assemblage (Consistency)
        for(int alpha=0;alpha<2;alpha++){
          for(int beta=0;beta<2;beta++){
            double na_mean=0.5*(nhat_p[alpha](beta)+nhat_m[alpha](beta));
            for(int k=0;k<3;k++){
              for(int j=0;j<nbFF_m;j++)
                m(j+k*nbFF_m,0)+=-(na_mean*Val_m[j]*(lbs->getphi0(beta,k))*(-lbs->getphi0(1,alpha)));
              for(int j=0;j<nbFF_p;j++)
                m(j+k*nbFF_p+nbdof_m,0)+=(na_mean*Val_p[j]*(lbs->getphi0(beta,k))*(-lbs->getphi0(1,alpha)));
            }
          }
        }
        // compute jump of u
        displacementjump(Val_m,nbFF_m,Val_p,nbFF_p,disp,ujump);
        Hhat_m->hat(lb_m,Cnt,nu);
        Hhat_p->hat(lb_p,Cnt,nu);
        Hhatmean->mean(Hhat_m,Hhat_p);
        // Assemblage (Compatibility)
        for(int alpha=0;alpha<2;alpha++){
          for(int j=0;j<nbFF_m;j++)
            for(int beta=0;beta<2;beta++)
              for(int gamma=0;gamma<2;gamma++)
                for(int delta=0;delta<2;delta++){
                  compC0PlateForceMembraneTerms(beta,gamma,delta,lb_m,lbs,Grads_m[j],ujump,me_comp);
                  for(int jj=0;jj<3;jj++)
                    m(j+jj*nbFF_m,0)+=(0.5*Hhat_m->get(alpha,beta,gamma,delta)*me_comp[jj]*(-lbs->getphi0(1,alpha)));
                }
          for(int j=0;j<nbFF_p;j++)
            for(int beta=0;beta<2;beta++)
              for(int gamma=0;gamma<2;gamma++)
                for(int delta=0;delta<2;delta++){
                  compC0PlateForceMembraneTerms(beta,gamma,delta,lb_p,lbs,Grads_p[j],ujump,me_comp);
                  for(int jj=0;jj<3;jj++)
                    m(j+jj*nbFF_p+nbdof_m,0)+=(0.5*Hhat_p->get(alpha,beta,gamma,delta)*me_comp[jj]*(-lbs->getphi0(1,alpha)));
                }
        }
        // Assemblage stability
        for(int alpha=0;alpha<2;alpha++)
            for(int beta=0;beta<2;beta++)
              for(int gamma=0;gamma<2;gamma++)
                for(int delta=0;delta<2;delta++){
                  stabC0PlateForceMembraneTerms(beta,delta,lbs,ujump,me_stab);
                  double temp = Hhatmean->get(alpha,beta,gamma,delta)*B2hs*(-lbs->getphi0(1,alpha))*(-lbs->getphi0(1,gamma));
                  for(int jj=0;jj<3;jj++){
                    for(int j=0;j<nbFF_m;j++)
                      m(j+jj*nbFF_m,0)+=-temp*Val_m[j]*me_stab[jj];
                    for(int j=0;j<nbFF_p;j++)
                      m(j+jj*nbFF_p+nbdof_m,0)+=temp*Val_p[j]*me_stab[jj];
                  }
                }

        // Out of Plane term
        Cst = weight*lbs->getJacobian()*Cs*B3hs; // Hooke tensor in "shearing" 1 component = Cs

        // compute B vector
        compute_Bs(lbs,Val_m,nbFF_m,Bs_m);
        compute_Bs(lbs,Val_p,nbFF_p,Bs_p);
        for(int j=0;j<nbFF_m;j++){
          stabilityC0PlateForceShearingTerms(Bs_m[j], Bs_m, Bs_p,nbFF_m,nbFF_p,disp,me_stab);
          for(int jj=0;jj<3;jj++)
            {m(j+jj*nbFF_m,0) += -Cst * me_stab[jj]; }//printf("-- %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);}
        }
        for(int j=0;j<nbFF_p;j++){
          stabilityC0PlateForceShearingTerms(Bs_p[j],Bs_m,Bs_p,nbFF_m,nbFF_p,disp,me_stab);
          for(int jj=0;jj<3;jj++)
            {m(j+(jj*nbFF_p)+nbdof_m,0) += Cst * me_stab[jj]; }//printf("++ %d %d %f\n",jj,kk,Cst * me_stab[jj][kk]);}

        }
        Val_m.clear(); Val_p.clear();
      }
    // Because component are push_back in Grads in gradfuvw idem for hess
    Grads_m.clear(); Grads_p.clear(); Hess_m.clear(); Hess_p.clear(); Grads.clear();
    }
  }
  else
    printf("not implemented\n");
}

void IsotropicElasticVirtualInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP,fullMatrix<double> &m)
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
    double Bhat_m[256][3][2][2];
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
      lb_m->set(velem[0],Grads_m); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lbs->set(iele,Grads,lb_m->gett0()); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??

      // PushForwardTensor
      lb_m->set_pushForward(lbs);

      /*printf("Local basis interface\n");
      lbs->print();
      printf("Local basis minus\n");
      lb_m->print();*/

      // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
      Compute_Bhat(lb_m,Hess_m,nbFF_m,Bhat_m);
      // Compute of Hooke hat tensor on minus and plus element
      Cmt = weight * lbs->getJacobian()  * Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
      Hhat_m->hat(lb_m,Cmt,nu);

      // Compute of Deltat_tilde
      compute_Deltat_tildeBound(lb_m,Grads_m,nbFF_m,Deltat_m,lbs);
      // loop on dof ATTENTION SAME NUMBER OF DOF for the two elements TODO take into account a different dof numbers between the two elements There ok because sym ??
      for(int j=0;j<nbFF_m;j++){
        for(int k=0;k<nbFF_m;k++){
          consC0PlateStiffnessTerms(Hhat_m,Bhat_m[j],Deltat_m[k],lbs,me_cons);
          compC0PlateStiffnessTerms(Hhat_m,Bhat_m[k],Deltat_m[j],lbs,me_comp);
          stabilityC0PlateStiffnessTerms(Hhat_m,Deltat_m[j], Deltat_m[k],lbs,me_stab);
          for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
              m(j+jj*nbFF_m,k+kk*nbFF_m) += -(me_cons[jj][kk] + me_comp[jj][kk] - Bhs * me_stab[jj][kk]);
        }
      }

      // Because component are push_back in Grads in gradfuvw idem for hess
      Grads_m.clear(); Hess_m.clear(); Grads.clear();
  }
//  m.print("InterfaceBound");
/*  // Compute the matrix by numerical perturbation Verif OK
  double eps=1.e-8;
  fullMatrix<double> fm(nbdof_m,1), fp(nbdof_m,1);
  fullMatrix<double> m2(nbdof_m,nbdof_m); m2.setAll(0.);
  fullMatrix<double> disp(nbdof_m,1);
  for(int j=0;j<nbdof_m;j++){
    disp.setAll(0.); fm.setAll(0.); fp.setAll(0.);
    disp(j,0)=eps;
    this->getInterForceOnBoundary(iele,npts,GP,disp,fp);
    disp(j,0)=-eps;
    this->getInterForceOnBoundary(iele,npts,GP,disp,fm);
    for(int k=0;k<nbdof_m;k++)
      m2(k,j)=(fp(k,0)-fm(k,0))/(2.*eps);
  }
  m2.print("Matrix by perturbation");*/
  }
  else
    printf("not implemented\n");
}

void IsotropicElasticVirtualInterfaceTermC0Plate::get(MElement *ele,int npts,IntPt *GP,const fullMatrix<double> &disp, fullMatrix<double> &m)
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
    m.resize(nbdof_m, 1);
    m.setAll(0.);
    double uem,vem,uep,vep;
    std::vector<TensorialTraits<double>::GradType> Grads_m;
    std::vector<TensorialTraits<double>::GradType> Grads;
    std::vector<TensorialTraits<double>::HessType> Hess_m;
    LocalBasis LBS,LBM;
    LocalBasis *lbs=&LBS; LocalBasis *lb_m=&LBM;
    double Bhat_m[256][3][2][2];
    LinearElasticShellHookeTensor HOOKehat_m;
    LinearElasticShellHookeTensor HOOKe_m;
    LinearElasticShellHookeTensor *Hhat_m=&HOOKehat_m;
    double Deltat_m[256][3][3];
    double Cmt;
    double me_cons[3];
    double me_comp[3];
    double me_stab[3];
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
      // Compute of gradient and hessian of shape functions on interface element and on elements
      // ATTENTION after multiplication by multipliers (functionspace 276) components are in the "second line of tensor change this ??
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(iele,u, v, w, Grads); // grad of shape fonction on interface element
      DgC0BilinearTerm<SVector3,SVector3>::space1.gradfuvw(velem[0],uem, vem, w, Grads_m); // w = 0
      DgC0BilinearTerm<SVector3,SVector3>::space1.hessfuvw(velem[0],uem, vem, w, Hess_m);

      // basis of elements and interface element
      lb_m->set(velem[0],Grads_m); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??
      lbs->set(iele,Grads,lb_m->gett0()); // This function can become a method of MElement Plante lors du calcul de l'énergie  à cause de  std::vector<TensorialTraits<SVector3>::GradType> prob de template??

      // PushForwardTensor
      lb_m->set_pushForward(lbs);

      // Compute of Bhat vector (1 component for now because 1 dof (z) ) --> it's a vector of length == nbFF
      Compute_Bhat(lb_m,Hess_m,nbFF_m,Bhat_m);
      // Compute of Hooke hat tensor on minus and plus element
      Cmt = weight * lbs->getJacobian()  * Cm; // Eh^3/(12(1-nu^2)) * weight gauss * jacobian
      Hhat_m->hat(lb_m,Cmt,nu);

      // Compute of Deltat_tilde
      compute_Deltat_tildeBound(lb_m,Grads_m,nbFF_m,Deltat_m,lbs);

      for(int j=0;j<nbFF_m;j++){
        consC0PlateForceTermsBound(Hhat_m,nbFF_m,Bhat_m[j],Deltat_m,lbs,disp,me_cons);
        compC0PlateForceTermsBound(nbFF_m,Hhat_m,Bhat_m,Deltat_m[j],lbs,disp,me_comp);
        stabilityC0PlateForceTermsBound(nbFF_m,Hhat_m,Deltat_m[j],Deltat_m,lbs,disp,me_stab);
        for(int jj=0;jj<3;jj++)
          m(j+jj*nbFF_m,0) += (me_cons[jj] - me_comp[jj] - Bhs * me_stab[jj] );
      }
      // Because component are push_back in Grads in gradfuvw idem for hess
      Grads_m.clear(); Hess_m.clear(); Grads.clear();
    }
  }
  else
    printf("not implemented\n");
};
