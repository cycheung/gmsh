// Gmsh - Copyright (C) 1997-2009 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

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

#if defined(HAVE_POST)
#include "PView.h"
#include "PViewData.h"
#endif

// DgC0PlateSolver
void DgC0PlateSolver::setMesh(const std::string meshFileName)
{
  pModel = new GModelWithInterface();
  pModel->readMSH(meshFileName.c_str());
  _dim = pModel->getNumRegions() ? 3 : 2;
  if (LagSpace) delete LagSpace;
  LagSpace=new DgC0LagrangeFunctionSpace(_tag);
}

void DgC0PlateSolver::readInputFile(const std::string fn)
{
  FILE *f = fopen(fn.c_str(), "r");
  char what[256];
  while(!feof(f)){
    if(fscanf(f, "%s", what) != 1) return;
    if (!strcmp(what, "ElasticPlaneStressLaw")){
      int num;
      double E,nu;
      if(fscanf(f, "%d %lf %lf",&num,&E,&nu) !=3) return;
      linearElasticLawPlaneStress *mlaw = new linearElasticLawPlaneStress(num,E,nu);
      this->addMaterialLaw(mlaw);
    }
    else if(!strcmp(what,"ElasticPlaneStressLawWithFracture")){
      int num;
      double E,nu,Gc,sigmac,beta,mu;
      if(fscanf(f, "%d %lf %lf %lf %lf %lf %lf",&num,&E,&nu,&Gc,&sigmac,&beta,&mu) !=7) return;
      linearElasticLawPlaneStressWithFracture *mlaw = new linearElasticLawPlaneStressWithFracture(num,E,nu,Gc,sigmac,beta,mu);
      this->addMaterialLaw(mlaw);
    }
    else if (!strcmp(what, "ElasticDomain")){
      dgLinearShellDomain* field = new dgLinearShellDomain(); // TODO creation of constructor
      int b,c,phys,matnum;
      double h;
      if(fscanf(f, "%d %d %lf %lf %d %d", phys, matnum, &b, &c) != 5) return;
      field->setTag(_tag);
      field->setFormulation(b); // Formulation of element cG/dG full dG
      field->setThickness(h);
      field->setPhysical(phys);
      field->setLawNum(matnum);
      // set the number of simpsons point for numerical integration on thickness (add one to the number if it is even)
      int msimp;
      c%2==0 ? msimp=++c : msimp=c;
      field->setNumberOfSimpsonPoints(msimp);
      field->g = new groupOfElements (_dim, phys);
      // push field in vector
      //elasticFields.push_back(field);
      domainVector.push_back(field);
    }
    else if (!strcmp(what, "Void")){}
    else if (!strcmp(what, "NodeDisplacement")){
      double val;
      int node, comp;
      if(fscanf(f, "%d %d %lf", &node, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (0, node);
      diri._f= simpleFunctionTime<double>(val);
      diri._comp=comp;
      diri._tag=node;
      diri.onWhat=BoundaryCondition::ON_VERTEX;
      allDirichlet.push_back(diri);
    }
    else if (!strcmp(what, "IndependentNodeDisplacement")){
      double val;
      int node, comp;
      if(fscanf(f, "%d %d %lf", &node, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (0, node);
      diri._f= simpleFunctionTime<double>(val,false);
      diri._comp=comp;
      diri._tag=node;
      diri.onWhat=BoundaryCondition::ON_VERTEX;
      allDirichlet.push_back(diri);
    }
    else if (!strcmp(what, "EdgeDisplacement")){
      double val;
      int edge, comp;
      if(fscanf(f, "%d %d %lf", &edge, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (1, edge);
      diri._f= simpleFunctionTime<double>(val);
      diri._comp=comp;
      diri._tag=edge;
      diri.onWhat=BoundaryCondition::ON_EDGE;
      allDirichlet.push_back(diri);
    }
    else if (!strcmp(what, "IndependentEdgeDisplacement")){
      double val;
      int edge, comp;
      if(fscanf(f, "%d %d %lf", &edge, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (1, edge);
      diri._f= simpleFunctionTime<double>(val,false);
      diri._comp=comp;
      diri._tag=edge;
      diri.onWhat=BoundaryCondition::ON_EDGE;
      allDirichlet.push_back(diri);
    }
    else if (!strcmp(what, "FaceDisplacement")){
      double val;
      int face, comp;
      if(fscanf(f, "%d %d %lf", &face, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (2, face);
      diri._f= simpleFunctionTime<double>(val);
      diri._comp=comp;
      diri._tag=face;
      diri.onWhat=BoundaryCondition::ON_FACE;
      allDirichlet.push_back(diri);
    }
    else if (!strcmp(what, "IndependentFaceDisplacement")){
      double val;
      int face, comp;
      if(fscanf(f, "%d %d %lf", &face, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (2, face);
      diri._f= simpleFunctionTime<double>(val,false);
      diri._comp=comp;
      diri._tag=face;
      diri.onWhat=BoundaryCondition::ON_FACE;
      allDirichlet.push_back(diri);
    }
    else if (!strcmp(what, "NodeForce")){
      double val1, val2, val3;
      int node;
      if(fscanf(f, "%d %lf %lf %lf", &node, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (0, node);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3));
      neu._tag=node;
      neu.onWhat=BoundaryCondition::ON_VERTEX;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "IndependentNodeForce")){
      double val1, val2, val3;
      int node;
      if(fscanf(f, "%d %lf %lf %lf", &node, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (0, node);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3),false);
      neu._tag=node;
      neu.onWhat=BoundaryCondition::ON_VERTEX;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "EdgeForce")){
      double val1, val2, val3;
      int edge;
      if(fscanf(f, "%d %lf %lf %lf", &edge, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (1, edge);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3));
      neu._tag=edge;
      neu.onWhat=BoundaryCondition::ON_EDGE;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "IndependentEdgeForce")){
      double val1, val2, val3;
      int edge;
      if(fscanf(f, "%d %lf %lf %lf", &edge, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (1, edge);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3),false);
      neu._tag=edge;
      neu.onWhat=BoundaryCondition::ON_EDGE;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "FaceForce")){
      double val1, val2, val3;
      int face;
      if(fscanf(f, "%d %lf %lf %lf", &face, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (2, face);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3));
      neu._tag=face;
      neu.onWhat=BoundaryCondition::ON_FACE;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "IndependentFaceForce")){
      double val1, val2, val3;
      int face;
      if(fscanf(f, "%d %lf %lf %lf", &face, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (2, face);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3),false);
      neu._tag=face;
      neu.onWhat=BoundaryCondition::ON_FACE;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "VolumeForce")){
      double val1, val2, val3;
      int volume;
      if(fscanf(f, "%d %lf %lf %lf", &volume, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (3, volume);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3));
      neu._tag=volume;
      neu.onWhat=BoundaryCondition::ON_VOLUME;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "IndependentVolumeForce")){
      double val1, val2, val3;
      int volume;
      if(fscanf(f, "%d %lf %lf %lf", &volume, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (3, volume);
      neu._f= simpleFunctionTime<SVector3>(SVector3(val1, val2, val3),false);
      neu._tag=volume;
      neu.onWhat=BoundaryCondition::ON_VOLUME;
      neu.setGaussIntegrationRule();
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "MeshFile")){
      char name[245];
      if(fscanf(f, "%s", name) != 1) return;
      setMesh(name);
    }
    else if (!strcmp(what, "Theta")){
      int num, num_phys=0;
      // First component number of physical groups where theta is imposed (to 0 for now)
      fscanf(f,"%d",&num);
      for(int i=0;i<num;i++){
        fscanf(f,"%d",&num_phys);
        std::vector<MVertex*> vv;
        pModel->getMeshVerticesForPhysicalGroup(1,num_phys,vv);
        pModel->insertThetaBound(num_phys,vv);
      }
    }
    else if(!strcmp(what, "Broken")){
      int numphys;
      fscanf(f,"%d",&numphys);
      this->addPhysInitBroken(numphys);
    }
    else if(!strcmp(what, "Solver")){
      int sol,sch;
      fscanf(f,"%d %d",&sol,&sch);
      this->setSolver(sol);
      this->setScheme(sch);
    }
    else if(!strcmp(what, "StabilityParameters")){
      double b1,b2,b3;
      fscanf(f,"%lf %lf %lf",&b1,&b2,&b3);
      this->setStabilityParameters(b1,b2,b3);
    }
    else if(!strcmp(what, "StaticNonLinearData")){
      int ns;
      double et,tol;
      fscanf(f,"%d %lf %lf",&ns,&et,&tol);
      this->setSNLData(ns,et,tol);
    }
    else if(!strcmp(what, "Archive")){
      int na;
      fscanf(f,"%d",&na);
      this->setStepBetweenArchiving(na);
    }
    else if(!strcmp(what, "ArchivingEdgeForce")){
      int numphys,comp;
      fscanf(f,"%d %d",&numphys,&comp);
      this->addArchivingEdgeForce(numphys,comp);
    }
    else if(!strcmp(what, "ArchivingNodeDisplacement")){
      int num,comp;
      fscanf(f,"%d %d",&num,&comp);
      this->addArchivingNodeDisplacement(num,comp);
    }
    else {
      Msg::Error("Invalid input : %s", what);
//      return;
    }
  }
  fclose(f);
}

void DgC0PlateSolver::init(){
  // create interfaceElement
  this->createInterfaceElement();

  for(std::vector<partDomain*>::iterator it = domainVector.begin(); it!=domainVector.end(); ++it){
    partDomain *dom = *it;
    materialLaw* mlaw = this->getMaterialLaw(dom->getLawNum());
    dom->setMaterialLaw(mlaw);
    dom->setGaussIntegrationRule();
    if(dom->IsInterfaceTerms()) hasInterface=true;
  }
}

void DgC0PlateSolver::addDgLinearElasticShellDomain(dgLinearShellDomain* ED, const int f, const int d){
  ED->g = new groupOfElements (d, f);
  dgLinearShellDomain *field = new dgLinearShellDomain(*ED);
  hasInterface = true;
  domainVector.push_back(field);
}

void DgC0PlateSolver::addMaterialLaw(materialLaw* mlaw){
  maplaw.insert(std::pair<int,materialLaw*>(mlaw->getNum(),mlaw));
}

void DgC0PlateSolver::addLinearElasticLawPlaneStress(linearElasticLawPlaneStress* mlaw){
  this->addMaterialLaw(mlaw);
}
void DgC0PlateSolver::addLinearElasticLawPlaneStressWithFracture(linearElasticLawPlaneStressWithFracture *mlaw){
  this->addMaterialLaw(mlaw);
}

materialLaw* DgC0PlateSolver::getMaterialLaw(const int num){
  return (maplaw.find(num)->second);
}

void DgC0PlateSolver::addTheta(const int numphys){
  std::vector<MVertex*> vv;
  pModel->getMeshVerticesForPhysicalGroup(1,numphys,vv);
  pModel->insertThetaBound(numphys,vv);
}

void DgC0PlateSolver::addPhysInitBroken(const int numphys){
  initbrokeninter.push_back(numphys);
}

void DgC0PlateSolver::createInterfaceElement(){
  // creation only if there has interfaceterm
  if(hasInterface){
    // Compute and add interface element to the model
    // Loop on mesh element and store each edge (which will be distributed between InterfaceElement, Boundary InterfaceElement and VirtualInterfaceElement)
    // A tag including the vertex number is used to identified the edge common to two elements. In this case a interface element is created
    std::map<unsigned long int,Iedge> map_edge;
    // loop on element
  //  for (int i = 0; i < elasticFields.size(); ++i)
    for(std::vector<partDomain*>::iterator itdom=domainVector.begin(); itdom !=domainVector.end(); ++itdom)
    {
      dgPartDomain* dom = dynamic_cast<dgPartDomain*>(*itdom);
      for (groupOfElements::elementContainer::const_iterator it = dom->g->begin(); it != dom->g->end(); ++it)
      {
        MElement *e=*it;
        int nedge = e->getNumEdges();
        for(int j=0;j<nedge;j++){
          std::vector<MVertex*> vv;
          e->getEdgeVertices(j,vv);
          Iedge ie = Iedge(vv,e,dom->getPhysical());
          unsigned long int key = ie.getkey();
          const std::map<unsigned long int,Iedge>::iterator it_edge=map_edge.find(key);
          if(it_edge == map_edge.end()) // The edge doesn't exist -->inserted into the map
            map_edge.insert(std::pair<unsigned long int,Iedge>(key,ie));
          else{ // create an interface element and remove the entry in map_edge
            MInterfaceElement interel =  MInterfaceElement(vv, 0, 0, e, it_edge->second.getElement()); // How to have num and part ?? TODO
            pModel->storeInterfaceElement(dom->getPhysical(),interel);
            map_edge.erase(it_edge);
          }
        }
      }
      pModel->getInterface(dom->getPhysical(),dom->gi); // Avoid this step ??
    }
    // At this stage the edge in map_edge must be separated into BoundaryInterfaceElement and Virtual InterfaceElement
    for(std::map<unsigned long int,Iedge>::iterator it = map_edge.begin(); it!=map_edge.end();++it){
      Iedge ie=it->second;
      std::vector<MVertex*> Mv= ie.getVertices();
      MInterfaceElement interel = MInterfaceElement(Mv,0,0,ie.getElement(),ie.getElement());
      pModel->storeVirtualInterfaceElement(interel);
    }
    std::vector<MInterfaceElement*> vie;
    pModel->getVirtualInterface(vie);

    std::vector<MVertex*> vv;
    std::vector<int> thetaBound;
    pModel->getphysBound(thetaBound);
    for(int i=0;i<thetaBound.size();i++){
      pModel->getThetaBound(thetaBound[i],vv);
      for(std::map<unsigned long int,Iedge>::iterator it_edge = map_edge.begin(); it_edge!=map_edge.end();++it_edge){
        for(int j=0;j<vv.size();j++){
          //loop on component of map_edge
          if(vv[j] == it_edge->second.getFirstInteriorVertex()){ // Ok because 2nd degree min BoundaryInterfaceElement is created only for this vertex
            std::vector<MVertex*> Mv = it_edge->second.getVertices();
            MInterfaceElement interel = MInterfaceElement(Mv, 0, 0, it_edge->second.getElement(),it_edge->second.getElement());
            pModel->storeBoundaryInterfaceElement(it_edge->second.getPhys(),interel);
            map_edge.erase(it_edge);
            break;
          }
        }
      }
      vv.clear();
    }
  //  for(int i=0;i<elasticFields.size();i++)
    for(std::vector<partDomain*>::iterator it = domainVector.begin(); it!=domainVector.end(); ++it){
      dgPartDomain *dgdom = dynamic_cast<dgPartDomain*>(*it);
      pModel->getBoundInterface(dgdom->getPhysical(),dgdom->gib);
    }
  }
}

void DgC0PlateSolver::createInterfaceElement_2(){
  // The contructor of dgGroupCollection create automatically the GroupsOfElements
  _groups = dgGroupCollection(this->pModel,this->_dim,1); // TODO ?? Add parameter to model to store order
  _groups.buildGroupsOfInterfaces();
  // Affichage temporaire pour vérification
/*  int nn = _groups.getNbFaceGroups();
  printf("Number of group of faces : %d\n",nn);
  for(int i=0;i<nn;i++){
    printf("Group of face number %d\n",i);
    dgGroupOfFaces *inter = _groups.getFaceGroup(i);
    int nnn = inter->getNbGroupOfConnections();
    printf("Number of connection group %d \n",nnn);
    for(int j=0;j<nnn;j++){
      const dgGroupOfConnections connec = inter->getGroupOfConnections(j);
      printf("Connection's group number %d\n",j);
      int nnnn = connec.getNbClosures();
      printf("Number of closures %d\n",nnnn);
      for(int k=0;k<nnnn;k++){
        printf("Closure number %d\n",k);
        std::vector<int> vec = connec.getClosure(k);
        for(int kk=0;kk<vec.size();kk++){
          printf(" %d ",vec[kk]);
        }
        printf("\n");
      }
    }
  }*/
}

void DgC0PlateSolver::solve()
{
// init data
this->init();
linearSystem<double> *lsys;
if(DgC0PlateSolver::whatSolver == Taucs){
    #if defined(HAVE_TAUCS)
      lsys = new linearSystemCSRTaucs<double>;
      printf("Taucs is chosen to solve\n");
    #else
      lsys = new linearSystemGmm<double>;
      lsys = dynamic_cast<linearSystemGmm<double>*>(lsys);
      dynamic_cast<linearSystemGmm<double>*>(lsys)->setNoisy(2);
      printf("Taucs is not installed\n Gmm is chosen to solve\n");
    #endif
}
else if(DgC0PlateSolver::whatSolver == Petsc){
    #if defined(HAVE_PETSC)
      lsys = new linearSystemPETSc<double>;
      printf("PETSc is chosen to solve\n");
    #else
      lsys = new linearSystemGmm<double>;
      lsys = dynamic_cast<linearSystemGmm<double>*>(lsys);
      dynamic_cast<linearSystemGmm<double>*>(lsys)->setNoisy(2);
      printf("PETSc is not installed\n Gmm is chosen to solve\n");
    #endif
}
else{
  lsys = new linearSystemGmm<double>;
  dynamic_cast<linearSystemGmm<double>*>(lsys)->setNoisy(2);
  printf("Gmm is chosen to solve\n");
}

  if (pAssembler) delete pAssembler;
  pAssembler = new dofManager<double>(lsys);

  // we first do all fixations. the behavior of the dofManager is to
  // give priority to fixations : when a dof is fixed, it cannot be
  // numbered afterwards

  // ATTENTION The BC must be rewrite to take into account different fields (CG/DG and fullDG field for exemple)
  std::cout <<  "Dirichlet BC"<< std::endl;
  std::vector<MInterfaceElement*> vinter;
  pModel->getVirtualInterface(vinter); // vector needed to impose boundary condition for fullDg formulation
  std::vector<MInterfaceElement*> vinternalInter;
  pModel->getInterface(vinternalInter);

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
          FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,vinter,vinternalInter);
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

  // displacement field
  displacementField ufield(pAssembler,domainVector,3,LagSpace->getId(),anoded);
  // Store stress and deformation at each gauss point
  // IPState "declaration" reserve place for data
  IPField<partDomain*,DgC0FunctionSpace<SVector3> > ipf(&domainVector,pAssembler,LagSpace, pModel, &ufield); // Todo fix this
  ipf.compute1state(IPState::initial);
  for(std::vector<partDomain*>::iterator itdom = domainVector.begin(); itdom!=domainVector.end(); ++itdom){
    partDomain *dom = *itdom;
    if(dom->getSolElemType() == SolElementType::ShellPlaneStress){
      dgLinearShellDomain *dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      dglshdom->InitializeTerms(*LagSpace,&ufield,&ipf,_beta1,_beta2,_beta3);
    }
    else if(dom->getSolElemType() == SolElementType::ShellPlaneStressWF){
      dgLinearShellDomain *dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      dglshdom->InitializeTerms(*LagSpace,&ufield,&ipf,_beta1,_beta2,_beta3);
    }
    else if(dom->getSolElemType() == SolElementType::ShellPlaneStressWTI){
      dgLinearShellDomain *dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      dglshdom->InitializeTerms(*LagSpace,&ufield,&ipf,_beta1,_beta2,_beta3);
    }
    else Msg::Error("Terms are not defined for domain %d",dom->getTag());
  }

  // Now we start the assembly process
  // First build the force vector
//  partDomain *dom0 = domainVector[0];
  if((dom0->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom0->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
     (dom0->getSolElemType()==SolElementType::ShellPlaneStressWF)){
    dgLinearShellDomain *dglsh0 = dynamic_cast<dgLinearShellDomain*>(dom0);
    // Integratioon Rule of ElasticField 0 for now fix it TODO ??
    std::cout <<  "Neumann BC"<< std::endl;
    for (unsigned int i = 0; i < allNeumann.size(); i++)
    {
      DgC0LoadTerm<SVector3> Lterm(*LagSpace,allNeumann[i]._f);
      if(!dglsh0->getFormulation())     // Use formulation of first field CHANGE THIS
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
    Msg::Error("External forces computation is not yet implemented for element without interface (dof numbering)");
  }

  // bulk material law
/*  for(std::vector<partDomain*>::iterator itdom = domainVector.begin(); itdom!=domainVector.end();++itdom)
  {
    // print the chosen formulation
    partDomain *dom = *itdom;
    // print for shell element
    if((dom->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
       (dom->getSolElemType()==SolElementType::ShellPlaneStressWF)){
      dgLinearShellDomain* dglshdom = dynamic_cast<dgLinearShellDomain*>(dom);
      if(dglshdom->getFormulation()) printf("Full Dg formulation is chosen for Domain %d\n",dglshdom->getTag());
      else printf("Cg/Dg formulation is chosen for Domain %d\n",dglshdom->getTag());
      printf("Number of integration point on thickness : %d\n",dglshdom->getmsimp());
    }

    // Initialization of elementary terms in function of the field and space
    //IsotropicElasticStiffBulkTermC0Plate Eterm(*LagSpace,elasticFields[i].getMaterialLaw(),elasticFields[i].getThickness(),elasticFields[i].getFormulation(),
    //                                           &ufield,&ipf,elasticFields[i].getSolElemType());
    // Assembling loop on Elementary terms
    MyAssemble(*(dom->getBilinearBulkTerm()),*LagSpace,dom->g->begin(),dom->g->end(),*(dom->getBulkGaussIntegrationRule()),
                 *pAssembler);

    if(dom->IsInterfaceTerms()){
      dgPartDomain* dgdom = dynamic_cast<dgPartDomain*>(dom);
      // Initialization of elementary  interface terms in function of the field and space
      //IsotropicElasticStiffInterfaceTermC0Plate IEterm(*LagSpace,elasticFields[i].getMaterialLaw(),_beta1,
      //                                            _beta2,_beta3,
      //                                            elasticFields[i].getThickness(),&ufield, &ipf, elasticFields[i].getSolElemType(),
      //                                            elasticFields[i].getFormulation());
      // Assembling loop on elementary interface terms
      MyAssemble(*(dgdom->getBilinearInterfaceTerm()),*LagSpace,dgdom->gi.begin(),dgdom->gi.end(),*(dgdom->getInterfaceGaussIntegrationRule()),
                        *pAssembler); // Use the same GaussQuadrature rule than on the boundary

      // Initialization of elementary  interface terms in function of the field and space
      //IsotropicElasticStiffVirtualInterfaceTermC0Plate VIEterm(*LagSpace,elasticFields[i].getMaterialLaw(),_beta1,
      //                                            _beta2,_beta3,
      //                                            elasticFields[i].getThickness(),&ufield,&ipf,elasticFields[i].getSolElemType(),true,elasticFields[i].getFormulation());
      // Assembling loop on elementary boundary interface terms
      MyAssemble(*(dgdom->getBilinearVirtualInterfaceTerm()),*LagSpace,dgdom->gib.begin(),dgdom->gib.end(),
                 *(dgdom->getInterfaceGaussIntegrationRule()),*pAssembler); // Use the same GaussQuadrature rule than on the boundary
    }
  }*/
  this->computeStiffMatrix(lsys,pAssembler,&ufield,&ipf);
  printf("-- done assembling!\n");
  lsys->systemSolve();
  printf("-- done solving!\n");

  // compute stress after solve
  ufield.update();
  ipf.compute1state(IPState::current);
  // save solution (for visualisation)
  ipf.buildView(domainVector,0.,1,"VonMises",-1,false);
  ufield.buildView(domainVector,0.,1,"displacement",-1,false);
}

void DgC0PlateSolver::setSNLData(const int ns, const double et, const double reltol){
  if(whatScheme==StaticNonLinear){
    numstep=ns;endtime=et;_tol=reltol;
  }
  else{
    Msg::Error("Impossible to set data for Static Non Linear scheme because another is chosen to solve the problem");
  }
}

void DgC0PlateSolver::solveLUA(){
  switch(whatScheme){
    case StaticLinear :
      this->solve();
      break;
    case StaticNonLinear :
      this->solveSNL();
      break;
  }
}

void DgC0PlateSolver::addDisp(std::string onwhat, const int numphys, const int comp, const double value){
  dirichletBC diri;
  const std::string node("Node");
  const std::string edge("Edge");
  const std::string face("Face");
  if(onwhat==node){
    diri.g = new groupOfElements (0, numphys);
    diri.onWhat=BoundaryCondition::ON_VERTEX;
  }
  else if(onwhat==edge){
    diri.g = new groupOfElements (1, numphys);
    diri.onWhat=BoundaryCondition::ON_EDGE;
  }
  else if(onwhat==face){
    diri.g = new groupOfElements (2, numphys);
    diri.onWhat=BoundaryCondition::ON_FACE;
  }
  else Msg::Error("Impossible to prescribe a displacement on a %s\n",onwhat.c_str());
  diri._f= simpleFunctionTime<double>(value);
  diri._comp=comp;
  diri._tag=numphys;
  allDirichlet.push_back(diri);
}

void DgC0PlateSolver::addIndepDisp(std::string onwhat, const int numphys, const int comp, const double value){
  dirichletBC diri;
  const std::string node("Node");
  const std::string edge("Edge");
  const std::string face("Face");
  if(onwhat==node){
    diri.g = new groupOfElements (0, numphys);
    diri.onWhat=BoundaryCondition::ON_VERTEX;
  }
  else if(onwhat==edge){
    diri.g = new groupOfElements (1, numphys);
    diri.onWhat=BoundaryCondition::ON_EDGE;
  }
  else if(onwhat==face){
    diri.g = new groupOfElements (2, numphys);
    diri.onWhat=BoundaryCondition::ON_FACE;
  }
  else Msg::Error("Impossible to prescribe a displacement on a %s\n",onwhat.c_str());
  diri._f= simpleFunctionTime<double>(value,false);
  diri._comp=comp;
  diri._tag=numphys;
  allDirichlet.push_back(diri);
}

void DgC0PlateSolver::addForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval){
  neumannBC neu;
  const std::string node("Node");
  const std::string edge("Edge");
  const std::string face("Face");
  const std::string volume("Volume");
  if(onwhat==node){
    neu.g = new groupOfElements (0, numphys);
    neu.onWhat=BoundaryCondition::ON_VERTEX;
  }
  else if(onwhat==edge){
    neu.g = new groupOfElements (1, numphys);
    neu.onWhat=BoundaryCondition::ON_EDGE;
  }
  else if(onwhat==face){
    neu.g = new groupOfElements (2, numphys);
    neu.onWhat=BoundaryCondition::ON_FACE;
  }
  else if(onwhat==volume){
    neu.g = new groupOfElements (3, numphys);
    neu.onWhat=BoundaryCondition::ON_VOLUME;
  }
  else  Msg::Error("Impossible to prescribe a force on a %s\n",onwhat.c_str());
  neu._f= simpleFunctionTime<SVector3>(SVector3(xval, yval, zval));
  neu._tag=numphys;
  neu.setGaussIntegrationRule();
  allNeumann.push_back(neu);
}

void DgC0PlateSolver::addIndepForce(std::string onwhat, const int numphys, const double xval, const double yval, const double zval){
  neumannBC neu;
  const std::string node("Node");
  const std::string edge("Edge");
  const std::string face("Face");
  const std::string volume("Volume");
  if(onwhat==node){
    neu.g = new groupOfElements (0, numphys);
    neu.onWhat=BoundaryCondition::ON_VERTEX;
  }
  else if(onwhat==edge){
    neu.g = new groupOfElements (1, numphys);
    neu.onWhat=BoundaryCondition::ON_EDGE;
  }
  else if(onwhat==face){
    neu.g = new groupOfElements (2, numphys);
    neu.onWhat=BoundaryCondition::ON_FACE;
  }
  else if(onwhat==volume){
    neu.g = new groupOfElements (3, numphys);
    neu.onWhat=BoundaryCondition::ON_VOLUME;
  }
  else  Msg::Error("Impossible to prescribe a force on a %s\n",onwhat.c_str());
  neu._f= simpleFunctionTime<SVector3>(SVector3(xval, yval, zval),false);
  neu._tag=numphys;
  neu.setGaussIntegrationRule();
  allNeumann.push_back(neu);
}


void DgC0PlateSolver::addArchivingEdgeForce(const int numphys, const int comp){
  // get the node of the edge
  std::vector<MVertex*> vv;
  pModel->getMeshVerticesForPhysicalGroup(1,numphys,vv);

  std::vector<Dof> vdof;// all dof include in the edge

  if(hasInterface){
    // shell element
    // Test on domain 0 TODO change this ??
    partDomain *dom0 = domainVector[0];
    if((dom0->getSolElemType()==SolElementType::ShellPlaneStress) or  (dom0->getSolElemType()==SolElementType::ShellPlaneStressWTI) or
        (dom0->getSolElemType()==SolElementType::ShellPlaneStressWF)){
      // build the dof
      dgLinearShellDomain *dglshdom0 = dynamic_cast<dgLinearShellDomain*>(dom0);
      if(!dglshdom0->getFormulation()){ //cG/dG case
        for(int i=0;i<vv.size();i++)
          vdof.push_back(Dof(vv[i]->getNum(), DgC0PlateDof::createTypeWithThreeInts(comp,LagSpace->getId())));
      }
      else{
        // find elements of the vertex
        for(std::vector<partDomain*>::iterator itdom=domainVector.begin(); itdom!=domainVector.end();++itdom){
            dgLinearShellDomain* dglshdom = dynamic_cast<dgLinearShellDomain*>(*itdom);
            for(groupOfElements::elementContainer::const_iterator it = dglshdom->g->begin(); it != dglshdom->g->end(); ++it){
              MElement *e = *it;
              for(int j=0;j<e->getNumVertices();j++){
                for(int k=0;k<vv.size();k++)
                  if(e->getVertex(j) == vv[k]){
                    vdof.push_back(Dof(e->getNum(),DgC0PlateDof::createTypeWithThreeInts(comp,LagSpace->getId(),j)));
                }
              }
            }
        }
      }
      // keys = 10*numphys + comp otherwise no way to archive different components
      int key = 10*numphys+comp;
      aef[key] = vdof;
      aefvalue[key] = 0.;

      // remove old file (linux only ??)
      std::ostringstream oss;
      oss << numphys;
      std::string s = oss.str();
      oss.str("");
      oss << comp;
      std::string s2 = oss.str();
      std::string rfname = "rm force"+s+"comp"+s2+".csv";
      system(rfname.c_str());
    }
    else{
      Msg::Error("Archiving force is only implemented for dg shell\n");
    }
  }
  else{
    Msg::Error("has no rule exist to make dof for formulation without interface it is impossible to archiving a force\n");
  }
}

void DgC0PlateSolver::addArchivingNodeDisplacement(const int num, const int comp){
  // no distinction between cG/dG and full Dg formulation. class Displacement Field manage it
  anoded.push_back(Dof(num,DgC0PlateDof::createTypeWithThreeInts(comp,LagSpace->getId())));
  // remove old file (linux only ??)
  std::ostringstream oss;
  oss << num;
  std::string s = oss.str();
  oss.str("");
  oss << comp;
  std::string s2 = oss.str();
  std::string rfname = "rm NodalDisplacement"+s+"comp"+s2+".csv";
  system(rfname.c_str());
}

std::vector<Dof> DgC0PlateSolver::getDofArchForce(){
    std::vector<Dof> vDof;
    for(std::map<int,std::vector<Dof> >::iterator it = aef.begin(); it!=aef.end(); ++it)
      for(int i=0;i<it->second.size();i++)
        vDof.push_back(Dof(it->second[i].getEntity(),it->second[i].getType()));
  return vDof;
  }

// Lua interaction
void DgC0PlateSolver::registerBindings(binding *b)
{
  classBinding *cb = b->addClass<DgC0PlateSolver>("DgC0PlateSolver");
  cb->setDescription("Solver class to solve plate problem with Cg/Dg or full Dg formulation");
  methodBinding *cm;
  cm = cb->addMethod("Input", &DgC0PlateSolver::readInputFile);
  cm->setArgNames("fn",NULL);
  cm->setDescription("Use this to give a txt file with data");
  cm = cb->addMethod("CreateInterfaceElement", &DgC0PlateSolver::createInterfaceElement);
  cm->setDescription("Create all interface elements");
  cm = cb->addMethod("solve", &DgC0PlateSolver::solveLUA);
  cm->setDescription("This function solve the linear elastic problem");
  cm = cb->setConstructor<DgC0PlateSolver,int>();
  cm->setArgNames("number",NULL);
  cm->setDescription("Class of C0DgPlateSolver to solve a linear elastic plate problem with Cg/Dg or full Dg formulation");
  cm = cb->addMethod("setScheme", &DgC0PlateSolver::setScheme);
  cm->setArgNames("number",NULL);
  cm->setDescription("Strategy to solve problem Static linear=0 Static Non linear=1");
  cm = cb->addMethod("SNLData", &DgC0PlateSolver::setSNLData);
  cm->setArgNames("ns","et","reltol",NULL);
  cm->setDescription("Set data for Static Non Linear Scheme. First argument=number of time step. Second argument = end time. Third= relative tolerance for iterattion default 1.e-6 ");
  cm = cb->addMethod("whichSolver", &DgC0PlateSolver::setSolver);
  cm->setArgNames("s",NULL);
  cm->setDescription("Select the solver for resolution Gmm=0 (default) Taucs=1 PETsc=2");
  cm = cb->addMethod("readmsh", &DgC0PlateSolver::setMesh);
  cm->setArgNames("meshFileName",NULL);
  cm->setDescription("Read a mesh file");
  cm = cb->addMethod("stepBetweenArchiving", &DgC0PlateSolver::setStepBetweenArchiving);
  cm->setArgNames("na",NULL);
  cm->setDescription("Set the number of step between 2 archiving");
  cm = cb->addMethod("addDgLinearElasticShellDomain", &DgC0PlateSolver::addDgLinearElasticShellDomain);
  cm->setArgNames("ED","f","d",NULL);
  cm->setDescription("Add a elastic domain class Second argument is the physical number of the domain. Third is the dimension");
  cm = cb->addMethod("AddLinearElasticLawPlaneStress", &DgC0PlateSolver::addLinearElasticLawPlaneStress);
  cm->setArgNames("mlaw",NULL);
  cm->setDescription("Add a material law");
  cm = cb->addMethod("AddLinearElasticLawPlaneStressWithFracture", &DgC0PlateSolver::addLinearElasticLawPlaneStressWithFracture);
  cm->setArgNames("mlaw",NULL);
  cm->setDescription("Add a material law");
  cm = cb->addMethod("AddThetaConstraint", &DgC0PlateSolver::addTheta);
  cm->setArgNames("numphys",NULL);
  cm->setDescription("Add a constrain on rotation (=0) for a physical group");
  cm = cb->addMethod("prescribedDisplacement", &DgC0PlateSolver::addDisp);
  cm->setArgNames("onwhat","numphys","comp","value",NULL);
  cm->setDescription("Add a prescribed displacement. First argument value (string) : Node, Edge, Face");
  cm = cb->addMethod("prescribedForce", &DgC0PlateSolver::addForce);
  cm->setArgNames("onwhat","numphys","xval","yval","zval",NULL);
  cm->setDescription("Add a prescribed force. First argument value (string) : Node, Edge, Face, Volume");
  cm = cb->addMethod("independentPrescribedDisplacement", &DgC0PlateSolver::addIndepDisp);
  cm->setArgNames("onwhat","numphys","comp","value",NULL);
  cm->setDescription("Add a prescribed displacement independent of time. First argument value (string) : Node, Edge, Face");
  cm = cb->addMethod("independentPrescribedForce", &DgC0PlateSolver::addIndepForce);
  cm->setArgNames("onwhat","numphys","xval","yval","zval",NULL);
  cm->setDescription("Add a prescribed force independent of time. First argument value (string) : Node, Edge, Face, Volume");
  cm = cb->addMethod("ArchivingEdgeForce", &DgC0PlateSolver::addArchivingEdgeForce);
  cm->setArgNames("numphys","comp",NULL);
  cm->setDescription("Archive an force on an edge. First argument is the physical number of the edge the second is the comp x=0 y=1 z=2");
  cm = cb->addMethod("ArchivingNodalDisplacement", &DgC0PlateSolver::addArchivingNodeDisplacement);
  cm->setArgNames("num","comp",NULL);
  cm->setDescription("Archiving a nodal displacement. First argument is the number of node the second is the comp x=0 y=1, z=2");
  cm = cb->addMethod("stabilityParameters", &DgC0PlateSolver::setStabilityParameters);
  cm->setArgNames("b1","b2","b3",NULL);
  cm->setDescription("Beta value for stabilization. The first one only for Cg/Dg formulation");
  cm = cb->addMethod("broken", &DgC0PlateSolver::addPhysInitBroken);
  cm->setArgNames("numphys",NULL);
  cm->setDescription("Initial broken of interface given by a physical number");
  // delete other functions ??
/*  cm = cb->addMethod("AddNodalDisplacement", &DgC0PlateSolver::addNodalDisp);
  cm->setArgNames("node","comp","value",NULL);
  cm->setDescription("Add a nodal displacement");
  cm = cb->addMethod("AddEdgeDisplacement", & DgC0PlateSolver::addEdgeDisp);
  cm->setArgNames("edge","comp","value",NULL);
  cm->setDescription("Add an edge displacement");
  cm = cb->addMethod("AddFaceDisplacement", &DgC0PlateSolver::addFaceDisp);
  cm->setArgNames("face","comp","value",NULL);
  cm->setDescription("Add face displacement");
  cm = cb->addMethod("AddNodeForce", &DgC0PlateSolver::addNodeForce);
  cm->setArgNames("node","xval","yval","zval",NULL)
  cm->setDescription("Add a nodal force");
  cm = cb->addMethod("AddEdgeForce", &DgC0PlateSolver::addEdgeForce);
  cm->setArgNames("edge","xval","yval","zval",NULL);
  cm->setDescription("Add an edge force");
  cm = cb->addMethod("AddFaceForce", &DgC0PlateSolver::addFaceForce);
  cm->setArgNames("face","xval","yval","zval",NULL);
  cm->setDescription("Add a face force");
  cm = cb->addMethod("AddVolumeForce", &DgC0PlateSolver::addVolumeForce);
  cm->setArgNames("volume","xval","yval","zval",NULL);
  cm->setDescription("Add a volume force");*/
}

