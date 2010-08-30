//
// C++ Interface: terms
//
// Description: Files with definition of function : DgC0PlateSolver::solveSNL()
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <string.h>
#include "GmshConfig.h"
#include "DgC0PlateSolver.h"
#include "linearSystemCSR.h"
#include "linearSystemPETSc.h"
#include "linearSystemGMM.h"
#include "Numeric.h"
#include "GModelWithInterface.h"
#include "DgC0PlateTerms.h"
#include "solverAlgorithms.h"
#include "quadratureRules.h"
#include "DgC0PlateSolverField.h"
#include "DgC0PlateAlgorithms.h"
#include "MPoint.h"
#include "IPState.h"
#include "IPField.h"
#include "SimpleFunctionTime.h"
#include "displacementField.h"
#include "dofManagerNonLinearSystem.h"

#if defined(HAVE_POST)
#include "PView.h"
#include "PViewData.h"
#endif

// Function to compute the initial norm
double DgC0PlateSolver::computeNorm0(linearSystem<double> *lsys, dofManager<double> *pAssembler,
                                     displacementField *ufield, IPField<partDomain*,DgC0FunctionSpace<SVector3> > *ipf,
                                     std::vector<MInterfaceElement*> &vinter){
  double normFext,normFint;
  // Set all component of RightHandSide to zero (only RightHandSide is used to compute norm
  lsys->zeroRightHandSide();

  // fext norm
  // compute ext forces
  // Test on dom0 change this TODO ??
  partDomain *dom0 = domainVector[0];
  if((dom0->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom0->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
     (dom0->getSolElemType()==SolElementType::ShellPlaneStressWF)){
    dgLinearShellDomain *dglshdom0 = dynamic_cast<dgLinearShellDomain*>(dom0);
    for (unsigned int i = 0; i < allNeumann.size(); i++)
    {
      DgC0LoadTerm<SVector3> Lterm(*LagSpace,allNeumann[i]._f);
      if(!dglshdom0->getFormulation())     // Use formulation of first field CHANGE THIS
        Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),*(allNeumann[i].integ),*pAssembler,false);
      else{ // The boundary condition on face are computed separately (because of research of interfaceElement linked to the BC)
        if(allNeumann[i].onWhat == BoundaryCondition::ON_FACE)
          Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),*(allNeumann[i].integ),*pAssembler,true);
        else
          Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),*(allNeumann[i].integ),*pAssembler,vinter);
      }
    }
  }
  else{
    Msg::Error("Computation of external force is not implemented for this solElementType (dof numbering)\n");
  }
  // norm
  normFext = lsys->normInfRightHandSide();

  // Internal forces
  lsys->zeroRightHandSide();
  // compute internal forces
  for(std::vector<partDomain*>::iterator itdom = domainVector.begin(); itdom!=domainVector.end(); ++itdom)
  {
    partDomain *dom = *itdom;
    // Assembling loop on Elementary terms
    MyAssemble(*(dom->getLinearBulkTerm()),*LagSpace,dom->g->begin(),dom->g->end(),*(dom->getBulkGaussIntegrationRule()),*pAssembler);

    if(dom->IsInterfaceTerms()){
      dgPartDomain* dgdom = dynamic_cast<dgPartDomain*>(dom);
      // Assembling loop on elementary interface terms
      MyAssemble(*(dgdom->getLinearInterfaceTerm()),*LagSpace,dgdom->gi.begin(),dgdom->gi.end(),*(dgdom->getInterfaceGaussIntegrationRule()),
                        *pAssembler); // Use the same GaussQuadrature rule than on the boundary
      // Assembling loop on elementary boundary interface terms
      MyAssemble(*(dgdom->getLinearVirtualInterfaceTerm()),*LagSpace,dgdom->gib.begin(),dgdom->gib.end(),*(dgdom->getInterfaceGaussIntegrationRule()),
                        *pAssembler); // Use the same GaussQuadrature rule than on the boundary
    }
  }
    //norm
    normFint = lsys->normInfRightHandSide();

    // archive
    return normFext+normFint;
}
// compute Fext-Fint
double DgC0PlateSolver::computeRightHandSide(linearSystem<double> *lsys, dofManager<double> *pAssembler,
                                    displacementField *ufield, IPField<partDomain*,DgC0FunctionSpace<SVector3> > *ipf,
                                    std::vector<MInterfaceElement*> &vinter){
  lsys->zeroRightHandSide();
  dofManagerNLS<double> *pA2 = dynamic_cast<dofManagerNLS<double>*>(pAssembler);
  pA2->clearRHSfixed();
  for(std::vector<partDomain*>::iterator itdom=domainVector.begin(); itdom!=domainVector.end();++itdom)
  {
    partDomain *dom = *itdom;
    dom->inverseLinearTermSign(); // Because rightHandSide = Fext - Fint
    MyAssemble(*(dom->getLinearBulkTerm()),*LagSpace,dom->g->begin(),dom->g->end(),*(dom->getBulkGaussIntegrationRule()),*pAssembler);

    if(dom->IsInterfaceTerms()){
      dgPartDomain *dgdom = dynamic_cast<dgPartDomain*>(dom);
      // Assembling loop on elementary interface terms
      MyAssemble(*(dgdom->getLinearInterfaceTerm()),*LagSpace,dgdom->gi.begin(),dgdom->gi.end(),*(dgdom->getInterfaceGaussIntegrationRule()),
                        *pAssembler); // Use the same GaussQuadrature rule than on the boundary
      // Assembling loop on elementary boundary interface terms
      MyAssemble(*(dgdom->getLinearVirtualInterfaceTerm()),*LagSpace,dgdom->gib.begin(),dgdom->gib.end(),
                 *(dgdom->getInterfaceGaussIntegrationRule()),*pAssembler); // Use the same GaussQuadrature rule than on the boundary
    }
    dom->inverseLinearTermSign(); // Otherwise Fint has not the right side
  }

  // save Fint component to archive ( write in file when Fint = Fext)
  pA2->getForces(aef,aefvalue);
  // compute ext forces
  partDomain *dom0 = domainVector[0];
  if((dom0->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom0->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
     (dom0->getSolElemType()==SolElementType::ShellPlaneStressWF)){
    dgLinearShellDomain *dglshdom0 = dynamic_cast<dgLinearShellDomain*>(dom0);
    for (unsigned int i = 0; i < allNeumann.size(); i++)
    {
      DgC0LoadTerm<SVector3> Lterm(*LagSpace,allNeumann[i]._f);
      if(!dglshdom0->getFormulation())     // Use formulation of first field CHANGE THIS
        Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),*(allNeumann[i].integ),*pAssembler,false);
      else{ // The boundary condition on face are computed separately (because of research of interfaceElement linked to the BC)
        if(allNeumann[i].onWhat == BoundaryCondition::ON_FACE)
          Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),*(allNeumann[i].integ),*pAssembler,true);
        else
          Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),*(allNeumann[i].integ),*pAssembler,vinter);
      }
    }
  }
  else{
    Msg::Error("Computation of external force is not implemented for this solElementType (dof numbering)\n");
  }
  return lsys->normInfRightHandSide();

}

void DgC0PlateSolver::computeStiffMatrix(linearSystem<double> *lsys, dofManager<double> *pAssembler,
                                    displacementField *ufield, IPField<partDomain*,DgC0FunctionSpace<SVector3> > *ipf){

  for(std::vector<partDomain*>::iterator itdom=domainVector.begin(); itdom!=domainVector.end(); ++itdom)
  {
    partDomain* dom = *itdom;
    // Assembling loop on Elementary terms
    MyAssemble(*(dom->getBilinearBulkTerm()),*LagSpace,dom->g->begin(),dom->g->end(),*(dom->getBulkGaussIntegrationRule()),*pAssembler);

    if(dom->IsInterfaceTerms()){
      dgPartDomain *dgdom = dynamic_cast<dgPartDomain*>(dom);
      // Assembling loop on elementary interface terms
      MyAssemble(*(dgdom->getBilinearInterfaceTerm()),*LagSpace,dgdom->gi.begin(),dgdom->gi.end(),
                 *(dgdom->getInterfaceGaussIntegrationRule()),*pAssembler); // Use the same GaussQuadrature rule than on the boundary
     // Assembling loop on elementary boundary interface terms
      MyAssemble(*(dgdom->getBilinearVirtualInterfaceTerm()),*LagSpace,dgdom->gib.begin(),dgdom->gib.end(),
                 *(dgdom->getInterfaceGaussIntegrationRule()),*pAssembler); // Use the same GaussQuadrature rule than on the boundary
    }
  }
}

// Newton Raphson scheme to solve one step
void DgC0PlateSolver::NewtonRaphson(linearSystem<double> *lsys, dofManager<double> *pAssembler,displacementField *ufield,
                                    IPField<partDomain*,DgC0FunctionSpace<SVector3> > *ipf,
                                    std::vector<MInterfaceElement*> &vinter){
  // compute ipvariable
  ipf->compute1state(IPState::current);

  // Compute the norm0 : norm(Fint) + norm(Fext) for relative convergence
  double norm0 = this->computeNorm0(lsys,pAssembler,ufield,ipf,vinter);

  // Compute Right Hand Side (Fext-Fint)
  double normFinf = this->computeRightHandSide(lsys,pAssembler,ufield,ipf,vinter);

  // loop until convergence
  int iter=0;
  double relnorm = normFinf/norm0;
  printf("iteration n° %d : residu : %lf\n",iter,relnorm);
  while(relnorm>_tol){
    iter++;
    // Compute Stiffness Matrix
    this->computeStiffMatrix(lsys,pAssembler,ufield,ipf);

    // Solve KDu = Fext-Fint
    lsys->systemSolve();

    // update displacement
    ufield->update();

    // update ipvariable
    ipf->compute1state(IPState::current);
    //double a = 100.;
    //ufield->buildView(pModel,elasticFields,a);
    // new forces
    normFinf = this->computeRightHandSide(lsys,pAssembler,ufield,ipf,vinter);

    // Check of convergence
    relnorm = normFinf/norm0;
    printf("iteration n° %d : residu : %lf\n",iter,relnorm);
    lsys->zeroMatrix(); // here otherwise if it has been in computeMatrix the first iteration is very low HOW ??
  }
}

// For now no initial displacement
void DgC0PlateSolver::solveSNL()
{
  Msg::Info("SNL Data : nstep =%d endtime=%f",this->getNumStep(),this->getEndTime());
  // init
  this->init();
  // Select the solver for linear system Ax=b (KDeltau = Fext-Fint)
  linearSystem<double> *lsys;
  if(DgC0PlateSolver::whatSolver == Taucs){
    #if defined(HAVE_TAUCS)
      lsys = new linearSystemCSRTaucs<double>;
      Msg::Info("Taucs is chosen to solve\n");
    #else
      lsys = new linearSystemGmm<double>;
      lsys = dynamic_cast<linearSystemGmm<double>*>(lsys);
      dynamic_cast<linearSystemGmm<double>*>(lsys)->setNoisy(2);
      Msg::Error("Taucs is not installed\n Gmm is chosen to solve\n");
    #endif
  }
  else if(DgC0PlateSolver::whatSolver == Petsc){
    #if defined(HAVE_PETSC)
      lsys = new linearSystemPETSc<double>;
      Msg::Info("PETSc is chosen to solve\n");
    #else
      lsys = new linearSystemGmm<double>;
      lsys = dynamic_cast<linearSystemGmm<double>*>(lsys);
      dynamic_cast<linearSystemGmm<double>*>(lsys)->setNoisy(2);
      Msg::Error("PETSc is not installed\n Gmm is chosen to solve\n");
    #endif
  }
  else{
    lsys = new linearSystemGmm<double>;
    dynamic_cast<linearSystemGmm<double>*>(lsys)->setNoisy(2);
    Msg::Info("Gmm is chosen to solve\n");
  }

  // Creation of dof manager. The Dof where the displacement is prescribed are fixe to zero
  // for initialization and numerotation of Dof (t=0 at initialization --> prescribed displacement =0)
  for(unsigned int i=0;i<allDirichlet.size();i++) allDirichlet[i]._f.setTime(0.);
  //for(unsigned int i=0;i<allNeumann.size();i++) allNeumann[i]._f.setTime(0.);

  if (pAssembler) delete pAssembler;
  pAssembler = new dofManagerNLS<double>(lsys,this->getDofArchForce(),0.);

  // we first do all fixations. the behavior of the dofManager is to
  // give priority to fixations : when a dof is fixed, it cannot be
  // numbered afterwards

  // ATTENTION The BC must be rewriten to take into account different fields (CG/DG and fullDG field for exemple)
  std::vector<MInterfaceElement*> vinter;
  pModel->getVirtualInterface(vinter); // vector needed to impose boundary condition for fullDg formulation
  std::vector<MInterfaceElement*> vinternalInterface;
  pModel->getInterface(vinternalInterface);

  // test on dom 0 change this TODO ??
  partDomain *dom0 = domainVector[0];
  if((dom0->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom0->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
     (dom0->getSolElemType()==SolElementType::ShellPlaneStressWF)){
   dgLinearShellDomain* dglshdom0 = dynamic_cast<dgLinearShellDomain*>(dom0);
   for(unsigned int i = 0; i < allDirichlet.size(); i++)
    {
      DgC0PlateFilterDofComponent filter(allDirichlet[i]._comp);
      if(!dglshdom0->getFormulation())
        FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,false);
      else{ // BC on face are computed separately
        if(allDirichlet[i].onWhat == BoundaryCondition::ON_FACE)
          FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,true);
        else
          FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,vinter,vinternalInterface);
      }
    }
    // we number the dofs : when a dof is numbered, it cannot be numbered
    // again with another number.
    for(std::vector<partDomain*>::iterator itdom = domainVector.begin(); itdom!=domainVector.end(); ++itdom)
    {
      partDomain *dom = *itdom;
      // Use formulation of first field CHANGE THIS
      NumberDofs(*LagSpace, dom->g->begin(), dom->g->end(),*pAssembler,dglshdom0->getFormulation());
    }
    // total number of unkowns to allocate the system
    double nunk = pAssembler->sizeOfR();
    // allocate system
    lsys->allocate(nunk);
  }
  else{
    Msg::Error("Numbering of dof is not implemented for this element type");
  }

  // iterative scheme (loop on timestep)
  double curtime = 0.;
  double dt = double(endtime)/double(numstep);
  displacementField ufield(pAssembler,domainVector,3,LagSpace->getId(),anoded); // 3 components by nodes User choice ??
  IPField<partDomain*,DgC0FunctionSpace<SVector3> > ipf(&domainVector,pAssembler,LagSpace,pModel,&ufield); // Field for GaussPoint
  ipf.compute1state(IPState::initial);
  ipf.initialBroken(pModel, initbrokeninter);
  ipf.copy(IPState::initial,IPState::previous); // if initial stress previous must be initialized before computation
  // Initialization of terms
  for(std::vector<partDomain*>::iterator itdom = domainVector.begin(); itdom!=domainVector.end(); ++itdom){
    partDomain *dom = *itdom;
    if(dom->getSolElemType()==SolElementType::ShellPlaneStress ){
      dgLinearShellDomain *dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      dglshdom->InitializeTerms(*LagSpace,&ufield,&ipf,_beta1,_beta2,_beta3);
    }
    else if(dom->getSolElemType() == SolElementType::ShellPlaneStressWF){
      dgLinearShellDomain *dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      dglshdom->InitializeTerms(*LagSpace,&ufield,&ipf,_beta1,_beta2,_beta3);
    }
    else if(dom->getSolElemType() == SolElementType::ShellPlaneStressWTI ){
      dgLinearShellDomain *dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      dglshdom->InitializeTerms(*LagSpace,&ufield,&ipf,_beta1,_beta2,_beta3);
    }
    else
      Msg::Error("Terms are not defined for domain %d",dom->getTag());
  }


  for(int ii=0;ii<numstep;ii++){
    curtime+=dt;
    printf("t= %lf on %lf\n",curtime,endtime);

  if((dom0->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom0->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
     (dom0->getSolElemType()==SolElementType::ShellPlaneStressWF)){
    // Fixation (prescribed displacement)
    dgLinearShellDomain* dglsh0 = dynamic_cast<dgLinearShellDomain*>(dom0);
    for(unsigned int i = 0; i < allDirichlet.size(); i++)
    {
      allDirichlet[i]._f.setTime(curtime);
      DgC0PlateFilterDofComponent filter(allDirichlet[i]._comp);
      if(!dglsh0->getFormulation())
        FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,false);
      else{ // BC on face are computed separately
        if(allDirichlet[i].onWhat == BoundaryCondition::ON_FACE)
          FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,true);
        else
          FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,vinter,vinternalInterface);
      }
    }
  }
  else{
    Msg::Error("Impossible to fix dof for this solElemType");
  }
    ufield.updateFixedDof();
    for (unsigned int i = 0; i < allNeumann.size(); i++) allNeumann[i]._f.setTime(curtime);

    // Solve one step by NR scheme
    NewtonRaphson(lsys,pAssembler,&ufield,&ipf,vinter);

    //archive must be rewritten after discussion with christophe
    if((ii+1)%nsba == 0){
      ufield.buildView(domainVector,curtime,ii+1,"displacement",-1,false);
      ipf.buildView(domainVector,curtime,ii+1,"VonMises",-1,false);
      ipf.buildView(domainVector,curtime,ii+1,"sigmaxx",0,false);
      ipf.buildView(domainVector,curtime,ii+1,"sigmayy",1,false);
      ipf.buildView(domainVector,curtime,ii+1,"tauxy",3,false);
    }
    // Archiving crack tip position (works only if there is no crack branching) no works anymore
//    ipf.archCrackTipPosition(curtime);

    // test of fracture "in a balanced configuration"
    ipf.evalFracture(IPState::current);

    // Archiving
    // Edge displacement
    ufield.archiving(curtime);
    // Edge force value;
    FILE *fp;
    for(std::map<int,double>::iterator it=aefvalue.begin(); it!=aefvalue.end(); ++it){
      std::ostringstream oss;
      oss << (it->first)/10;
      std::string s = oss.str();
      oss.str("");
      oss << (it->first)%10;
      std::string s2 = oss.str();
      std::string fname = "force"+s+"comp"+s2+".csv";
      fp = fopen(fname.c_str(),"a");
      fprintf(fp,"%lf;%lf\n",curtime,it->second);
      fclose(fp);
    }

    // write M(Dr) and N(Dun) curves
/*    int ngauss = 4;
    int numinter = 15;
    double M, dr,delta,N,du;
    for(int i=0;i<ngauss;i++){
      ipf.getData(numinter,i,elasticFields[0].getMaterialLaw(),N,M,du,dr,delta);
      std::ostringstream oss;
      oss << i;
      std::string s = oss.str();
      std::string s1("MDelta");
      std::string ss1("NDelta");
      std::string s2(".csv");
      std::string fileName  = s1+s+s2;
      std::string fN2 = ss1+s+s2;
      FILE *fp = fopen(fileName.c_str(), "a");
      fprintf(fp,"%lf;%lf;%lf;%lf\n",curtime,delta,dr,M);
      fclose(fp);
      fp = fopen(fN2.c_str(), "a");
      fprintf(fp,"%lf;%lf;%lf;%lf\n",curtime,delta,du,N);
      fclose(fp);
    }*/

    // next step for ipvariable
    ipf.nextStep();
  }
}

