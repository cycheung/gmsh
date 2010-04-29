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

#if defined(HAVE_POST)
#include "PView.h"
#include "PViewData.h"
#endif

void DgC0PlateSolver::setMesh(const std::string &meshFileName)
{
  pModel = new GModelWithInterface();
  pModel->readMSH(meshFileName.c_str());
  _dim = pModel->getNumRegions() ? 3 : 2;
  if (LagSpace) delete LagSpace;
  // Faudra quelque chose dans le fichier de données qui permettra de choisir
  //if (_dim==3) LagSpace=new VectorLagrangeFunctionSpace(_tag);
  //if (_dim==2) LagSpace=new VectorLagrangeFunctionSpace(_tag,VectorLagrangeFunctionSpace::VECTOR_X,VectorLagrangeFunctionSpace::VECTOR_Y);
  // Plaque un ddl par noeud
  //LagSpace=new VectorLagrangeFunctionSpace(_tag,VectorLagrangeFunctionSpace::VECTOR_Z);
  // Plate 3 dof per node
  LagSpace=new DgC0LagrangeFunctionSpace(_tag);

}

void DgC0PlateSolver::readInputFile(const std::string &fn)
{
  FILE *f = fopen(fn.c_str(), "r");
  char what[256];
  while(!feof(f)){
    if(fscanf(f, "%s", what) != 1) return;
    if (!strcmp(what, "ElasticDomain")){
      DGelasticField field;
      int physical, b;
      if(fscanf(f, "%d %lf %lf %lf %lf %lf %lf %d", &physical, &field._E, &field._nu, &field._beta1, &field._beta2, &field._beta3, &field._h, &b) != 8) return;
      field._tag = _tag;
      field._phys= physical; // needed to impose BC on face
      field.setFormulation(b);
      field.setMaterialLaw(materialLaw::linearElasticPlaneStress); // Elastic Domain --> linear elastic law
      field.g = new groupOfElements (_dim, physical);
      // Add the Interface Elements
      field.gi = pModel->getInterface(physical); // TODO integrate with field.g
      elasticFields.push_back(field);
    }
    else if (!strcmp(what, "Void")){
      //      elasticField field;
      //      create enrichment ...
      //      create the group ...
      //      assign a tag
      //      elasticFields.push_back(field);
    }
    else if (!strcmp(what, "NodeDisplacement")){
      double val;
      int node, comp;
      if(fscanf(f, "%d %d %lf", &node, &comp, &val) != 3) return;
      dirichletBC diri;
      diri.g = new groupOfElements (0, node);
      diri._f= simpleFunction<double>(val);
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
      diri._f= simpleFunction<double>(val);
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
      diri._f= simpleFunction<double>(val);
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
      neu._f= simpleFunction<SVector3>(SVector3(val1, val2, val3));
      neu._tag=node;
      neu.onWhat=BoundaryCondition::ON_VERTEX;
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "EdgeForce")){
      double val1, val2, val3;
      int edge;
      if(fscanf(f, "%d %lf %lf %lf", &edge, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (1, edge);
      neu._f= simpleFunction<SVector3>(SVector3(val1, val2, val3));
      neu._tag=edge;
      neu.onWhat=BoundaryCondition::ON_EDGE;
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "FaceForce")){
      double val1, val2, val3;
      int face;
      if(fscanf(f, "%d %lf %lf %lf", &face, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (2, face);
      neu._f= simpleFunction<SVector3>(SVector3(val1, val2, val3));
      neu._tag=face;
      neu.onWhat=BoundaryCondition::ON_FACE;
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "VolumeForce")){
      double val1, val2, val3;
      int volume;
      if(fscanf(f, "%d %lf %lf %lf", &volume, &val1, &val2, &val3) != 4) return;
      neumannBC neu;
      neu.g = new groupOfElements (3, volume);
      neu._f= simpleFunction<SVector3>(SVector3(val1, val2, val3));
      neu._tag=volume;
      neu.onWhat=BoundaryCondition::ON_VOLUME;
      allNeumann.push_back(neu);
    }
    else if (!strcmp(what, "MeshFile")){
      char name[245];
      if(fscanf(f, "%s", name) != 1) return;
      setMesh(name);
    }
    else if (!strcmp(what, "Theta")){
      int num, num_phys=0;
      std::vector<MVertex*> vv;
      // First component number of physical groups where theta is imposed (to 0 for now)
      fscanf(f,"%d",&num);
      for(int i=0;i<num;i++){
        fscanf(f,"%d",&num_phys);
        // Store the group (the interfaceElement will be created later)
        pModel->getMeshVerticesForPhysicalGroup(1,num_phys,vv);
        pModel->insertThetaBound(num_phys,vv);
      }
      // store the boundary interface element to the elastifFields
      for(int i=0;i<elasticFields.size();i++) elasticFields[i].gib=pModel->getBoundInterface(i);
    }
    else if(!strcmp(what, "Solver")){
      int s;
      fscanf(f,"%d",&s);
      this->setSolver(s);
    }
    else {
      Msg::Error("Invalid input : %s", what);
//      return;
    }
  }
  fclose(f);
}

void DgC0PlateSolver::createInterfaceElement(){
  // Compute and add interface element to the model
  // Loop on mesh element and store each edge (which will be distributed between InterfaceElement, Boundary InterfaceElement and VirtualInterfaceElement)
  // A tag including the vertex number is used to identified the edge common to two elements. In this case a interface element is created
  std::map<unsigned long int,Iedge> map_edge;
  // loop on element
  for (int i = 0; i < elasticFields.size(); ++i)
  {
    for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
    {
      MElement *e=*it;
      int nedge = e->getNumEdges();
      for(int j=0;j<nedge;j++){
        std::vector<MVertex*> vv;
        e->getEdgeVertices(j,vv);
        Iedge ie = Iedge(vv,e,i);
        unsigned long int key = ie.getkey();
        const std::map<unsigned long int,Iedge>::iterator it_edge=map_edge.find(key);
        if(it_edge == map_edge.end()) // The edge doesn't exist -->inserted into the map
          map_edge.insert(std::pair<unsigned long int,Iedge>(key,ie));
        else{ // create an interface element and remove the entry in map_edge
          MInterfaceElement interel =  MInterfaceElement(vv, 0, 0, e, it_edge->second.getElement()); // How to have num and part ?? TODO
          pModel->storeInterfaceElement(elasticFields[i]._phys,interel);
          map_edge.erase(it_edge);
        }
      }
    }
    elasticFields[i].gi=pModel->getInterface(elasticFields[i]._phys); // Avoid this step ??
  }
  // At this stage the edge in map_edge must be separated into BoundaryInterfaceElement and Virtual InterfaceElement
  //pModel->generateInterfaceElementsOnBoundary(num_phys,elasticFields); //TODO don't pass elastic field but I don't know how to access to element in GModel ??
  for(std::map<unsigned long int,Iedge>::iterator it = map_edge.begin(); it!=map_edge.end();++it){
    Iedge ie=it->second;
    std::vector<MVertex*> Mv= ie.getVertices();
    MInterfaceElement interel = MInterfaceElement(Mv,0,0,ie.getElement(),ie.getElement());
    pModel->storeVirtualInterfaceElement(interel);
  }
  std::vector<MInterfaceElement*> vie = pModel->getVirtualInterface();

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
  for(int i=0;i<elasticFields.size();i++)
    elasticFields[i].gib=pModel->getBoundInterface(i);
  //for(int j=0;j<elasticFields[0].gib.size();j++)
  // printf("elem - %d elem + %d\n",elasticFields[0].gib[j]->getElem(0)->getNum(),elasticFields[0].gib[j]->getElem(1)->getNum());
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
  std::vector<MInterfaceElement*> vinter = pModel->getVirtualInterface(); // vector needed to impose boundary condition for fullDg formulation
  for (unsigned int i = 0; i < allDirichlet.size(); i++)
  {
    DgC0PlateFilterDofComponent filter(allDirichlet[i]._comp);
    if(!elasticFields[0].getFormulation())
      FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,false);
    else{ // BC on face are computed separately
      if(allDirichlet[i].onWhat == BoundaryCondition::ON_FACE)
        FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,true);
      else
        FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter,vinter);
    }
 }

  // we number the dofs : when a dof is numbered, it cannot be numbered
  // again with another number.
  for (unsigned int i = 0; i < elasticFields.size(); ++i)
  {
    // Use formulation of first field CHANGE THIS
    NumberDofs(*LagSpace, elasticFields[i].g->begin(), elasticFields[i].g->end(),*pAssembler,elasticFields[0].getFormulation());
  }

  // Now we start the assembly process
  // First build the force vector

  GaussQuadrature Integ_Boundary(GaussQuadrature::Val);
  GaussQuadrature Integ_Bulk(GaussQuadrature::GradGrad);
  std::cout <<  "Neumann BC"<< std::endl;
  for (unsigned int i = 0; i < allNeumann.size(); i++)
  {
    DgC0LoadTerm<SVector3> Lterm(*LagSpace,allNeumann[i]._f);
    if(!elasticFields[0].getFormulation())     // Use formulation of first field CHANGE THIS
      Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),Integ_Boundary,*pAssembler,false);
    else{ // The boundary condition on face are computed separately (because of research of interfaceElement linked to the BC)
      if(allNeumann[i].onWhat == BoundaryCondition::ON_FACE)
        Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),Integ_Boundary,*pAssembler,true);
      else
        Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),Integ_Boundary,*pAssembler,vinter);
    }
  }

  // Store stress and deformation at each gauss point
  // IPState "declaration" reserve place for data
  // Just storage allocation no initialization in AIPS constructor
  AllIPState AIPS(pModel, elasticFields, Integ_Bulk, Integ_Boundary);
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> > ipf(&elasticFields,pAssembler,LagSpace,
                                                           &Integ_Bulk, &Integ_Boundary, &AIPS);
  ipf.compute1state(IPState::initial);
  // depends on ipvar make a switch or something else ??
  // loop on all component of AIPS to initialize the variable at gauss point
 //  AIPS.copy(IPState::initial,IPState::current);

  // bulk material law
  for (unsigned int i = 0; i < elasticFields.size(); i++)
  {
    // print the chosen formulation
    if(elasticFields[i].getFormulation()) printf("Full Dg formulation is chosen for elasticField %d\n",elasticFields[i]._tag);
    else printf("Cg/Dg formulation is chosen for elasticField %d\n",elasticFields[i]._tag);
    // Initialization of elementary terms in function of the field and space
    IsotropicElasticBulkTermC0Plate Eterm(*LagSpace,elasticFields[i]._E,elasticFields[i]._nu,elasticFields[i]._h);
    // Assembling loop on Elementary terms
    AssembleBulk(Eterm,*LagSpace,elasticFields[i].g->begin(),elasticFields[i].g->end(),Integ_Bulk,
                 *pAssembler, elasticFields[i].getFormulation());

    // Initialization of elementary  interface terms in function of the field and space
    IsotropicElasticInterfaceTermC0Plate IEterm(*LagSpace,elasticFields[i]._E,elasticFields[i]._nu,elasticFields[i]._beta1,
                                                elasticFields[i]._beta2,elasticFields[i]._beta3,
                                                elasticFields[i]._h,elasticFields[i].getFormulation());
    // Assembling loop on elementary interface terms
    AssembleInterface(IEterm,*LagSpace,elasticFields[i].g,elasticFields[i].gi,Integ_Boundary,
                      *pAssembler,elasticFields[i].getFormulation()); // Use the same GaussQuadrature rule than on the boundary

    // Initialization of elementary  interface terms in function of the field and space
    IsotropicElasticVirtualInterfaceTermC0Plate VIEterm(*LagSpace,elasticFields[i]._E,elasticFields[i]._nu,elasticFields[i]._beta1,
                                                elasticFields[i]._beta2,elasticFields[i]._beta3,
                                                elasticFields[i]._h,elasticFields[i].getFormulation());
    // Assembling loop on elementary boundary interface terms
    AssembleInterface(VIEterm,*LagSpace,elasticFields[i].g,elasticFields[i].gib,Integ_Boundary,
                      *pAssembler,elasticFields[i].getFormulation()); // Use the same GaussQuadrature rule than on the boundary

  }
printf("-- done assembling!\n");
lsys->systemSolve();
printf("-- done solving!\n");

// compute stress after solve
//AIPS.nextStep();
ipf.compute1state(IPState::current);
//AIPS.stressAndDefo(IPState::current,elasticFields,pAssembler,LagSpace,Integ_Bulk,Integ_Boundary);
// print result
/*IntPt *GP;
IPState::whichState ws = IPState::current;
      for(int i=0;i<elasticFields.size();i++){
        if(!elasticFields[i].getFormulation()){ //cg/dg formulation
          // The object is initialized thanks to interfaceElement. Other Possibilities ?? if CG can be done with element but not for cg/dg and full dg
          // loop
          for(std::vector<MInterfaceElement*>::iterator it=elasticFields[i].gi.begin(); it!=elasticFields[i].gi.end();++it){
            MInterfaceElement *ie = *it;
            // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
            MElement *em = ie->getElem(0);
            int edge = ie->getEdgeNumber(0);
            int npts_inter=Integ_Boundary.getIntPoints(ie,&GP);
            for(int j=0;j<npts_inter;j++){
              IPnum key=IPnum(em->getNum(),IPnum::createTypeWithTwoInts(edge,j));
              IPState *ips = AIPS.getIPstate(&key);
              IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
              printf("elem %d edge %d gaussPoint %d sxx =%f\n",em->getNum(),edge,j,ipv->getSigma(component::xx));
            }
          }
        }
        else{ //full dg formulation
          // loop
          for(std::vector<MInterfaceElement*>::iterator it=elasticFields[i].gi.begin(); it!=elasticFields[i].gi.end();++it){
            MInterfaceElement *ie = *it;
            // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
            MElement *em = ie->getElem(0);
            MElement *ep = ie->getElem(1);
            int edgem = ie->getEdgeNumber(0);
            int edgep = ie->getEdgeNumber(1);
            int npts_inter=Integ_Boundary.getIntPoints(ie,&GP);
            for(int j=0;j<npts_inter;j++){
              IPnum key=IPnum(em->getNum(),IPnum::createTypeWithTwoInts(edgem,j));
              IPState *ips = AIPS.getIPstate(&key);
              IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
              printf("elem %d edge %d gaussPoint %d sxx =%f\n",em->getNum(),edgem,j,ipv->getSigma(component::xx));
              IPnum keyp=IPnum(ep->getNum(),IPnum::createTypeWithTwoInts(edgep,j));
              ips = AIPS.getIPstate(&keyp);
              ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
              printf("elem %d edge %d gaussPoint %d sxx =%f\n",ep->getNum(),edgep,j,ipv->getSigma(component::xx));
            }
          }
        }
        for(std::vector<MInterfaceElement*>::iterator it=elasticFields[i].gib.begin(); it!=elasticFields[i].gib.end();++it){
          MInterfaceElement *ie = *it;
          // get info for element - and + (gauss point on interface are only created for element - has cg/dg)
          MElement *em = ie->getElem(0);
          int edge = ie->getEdgeNumber(0);
          int npts_inter=Integ_Boundary.getIntPoints(ie,&GP);
          for(int j=0;j<npts_inter;j++){
            IPnum key=IPnum(em->getNum(),IPnum::createTypeWithTwoInts(edge,j));
            IPState *ips = AIPS.getIPstate(&key);
            //ips->getState(IPState::initial);
            //IPVariablePlate *ipv = <dynamic_cast IPVariablePlate *> ips->getState(IPState::initial);
            IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
            printf("elem %d edge %d gaussPoint %d sxx =%f\n",em->getNum(),edge,j,ipv->getSigma(component::xx));
          }
        }

        for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it){
          MElement *e = *it;
          int edge = e->getNumEdges();
          int npts_bulk=Integ_Bulk.getIntPoints(e,&GP);
          for(int j=0;j<npts_bulk;j++){
            IPnum key=IPnum(e->getNum(),IPnum::createTypeWithTwoInts(edge,j));
            IPState *ips = AIPS.getIPstate(&key);
            IPVariablePlate *ipv = dynamic_cast<IPVariablePlate*>(ips->getState(ws));
            printf("elem %d edge %d gaussPoint %d sxx =%f\n",e->getNum(),edge,j,ipv->getSigma(component::xx));
          }
        }
      }*/
 // Save solution
/*FILE *fp1 = fopen("matrix.txt","w");
FILE *fp2 = fopen("force.txt","w");
FILE *fp3 = fopen("sol.txt","w");
for(int i=0;i<48;i++){
  for(int j=0;j<48;j++)
    fprintf(fp1,"%f ",lsys->getFromMatrix(i,j));
  fprintf(fp2, "%f\n",lsys->getFromRightHandSide(i));
  fprintf(fp3, "%f\n",lsys->getFromSolution(i));
  fprintf(fp1,"\n");
}
fclose(fp1); fclose(fp2); fclose(fp3);*/
// Plante Bug problème de template ??
/*  double energ=0;
  for (unsigned int i = 0; i < elasticFields.size(); i++)
  {
    SolverField<SVector3> Field(pAssembler, LagSpace);
    IsotropicElasticTermC0Plate Eterm(Field,elasticFields[i]._E,elasticFields[i]._nu,elasticFields[i]._h);
    BilinearTermToScalarTerm<SVector3,SVector3> Elastic_Energy_Term(Eterm);
    Assemble(Elastic_Energy_Term,elasticFields[i].g->begin(),elasticFields[i].g->end(),Integ_Bulk,energ);
  }
  printf("elastic energy=%f\n",energ);*/
}



#if defined(HAVE_POST)

// Function returning the type of element (Implemented only for supported element by FullDgPlate
int whichType(MElement *e){
  int type;
  int numEdge = e->getNumEdges();
  int order = e->getPolynomialOrder();
  switch(numEdge){
    case 2 :
      switch(order){
        case 1 :
          type=1; break;
        case 2 :
          type=8; break;
        case 3 :
          type=26; break;
        default : printf("Warning unknow element type in builDisplacementView"); type=0; break;
      } break;
    case 3 :
      switch(order){
        case 1 :
          type=2; break;
        case 2 :
          type=9; break;
        case 3 :
          type=20; break;
        default :  printf("Warning unknow element type in builDisplacementView"); type=0; break;
      } break;
    case 4 :
      switch(order){
        case 1 :
          type=3; break;
        case 2 :
          type=16; break;
        case 3 :
          type=39 ; break;
        default :  printf("Warning unknow element type in builDisplacementView"); type=0; break;
      } break;
    default : printf("Warning unknow element type in builDisplacementView"); type=0; break;

  }

  return type;
}

PView* DgC0PlateSolver::buildDisplacementView (const std::string &postFileName)
{
  if(!elasticFields[0].getFormulation()){ // No change if CG/DG change this because test on elasticFields[0] only
    std::set<MVertex*> v;
    for (unsigned int i = 0; i < elasticFields.size(); ++i)
    {
      for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
      {
        MElement *e=*it;
        for (int j = 0; j < e->getNumVertices(); ++j) v.insert(e->getVertex(j));
      }
    }
    std::map<int, std::vector<double> > data;
    DgC0SolverField<SVector3> Field(pAssembler, LagSpace);
    for ( std::set<MVertex*>::iterator it = v.begin(); it != v.end(); ++it)
    {
      SVector3 val;
      MPoint p(*it);
      Field.getVertexDisplacement(&p,val,false,0);
      std::vector<double> vec(3);vec[0]=val(0);vec[1]=val(1);vec[2]=val(2);
      data[(*it)->getNum()]=vec;
    }
    PView *pv = new PView (postFileName, "NodeData", pModel, data, 0.0);
    return pv;
  }
  else{ // Full Dg case all files is save here and function pv->getData->writeMSH is not used so return null pointer
    std::set<MElement*> v; //set with element has each element has got its own dofs
    for (unsigned int i = 0; i < elasticFields.size(); ++i)
      for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
        v.insert(*it);

    SVector3 val(0.,0.,0.);
    std::map<unsigned long int,MVertex*> map_vertex; // map with the number of node associated to a vertex
    std::map<MElement*,std::vector<unsigned long int> > map_elem;
    std::map<unsigned long int,SVector3> data;
    unsigned long int num=1;
    DgC0SolverField<SVector3> Field(pAssembler, LagSpace);
    std::vector<unsigned long int> num_vertex;
    for( std::set<MElement*>::iterator it = v.begin(); it != v.end() ;++it){
      MElement *e=*it;
      for(int j=0;j<e->getNumVertices();j++){
        Field.getVertexDisplacement(e,val,true,j);
        map_vertex[num]=e->getVertex(j);
        data[num]=val;
        num_vertex.push_back(num);
        num++;
      }
    //map_elem.insert(std::pair<int,std::vector<unsigned long int> >(e->getNum(),num_vertex));
    map_elem[e] = num_vertex;
    num_vertex.clear();
  }

  // Now all informations are disponible for write GMSH file make a function in PViewData writeMSH with this code (will be necessary to view stress or strain)
  bool binary=false;               // arg of function
  std::string fileName = "dispDg.msh"; // arg of function
  std::string name_t = "displacement"; //arg of function
  double time = 0.; // arg of function
  int ts =0; //arg of function What is ts;

  const int numComp=3; // not arg because SVector3
  FILE *fp = fopen(fileName.c_str(), "w");
  if(!fp){
    Msg::Error("Unable to open file '%s'", fileName.c_str());
    return false;
  }

  if(binary) Msg::Warning("Binary write not implemented yet");

  fprintf(fp, "$MeshFormat\n2 0 8\n$EndMeshFormat\n");
  fprintf(fp, "$Nodes\n");
  fprintf(fp, "%d\n", (int)map_vertex.size());
  for(std::map<unsigned long int, MVertex*>::iterator it = map_vertex.begin();it!=map_vertex.end();++it){
      fprintf(fp, "%d %.16g %.16g %.16g\n", (*it).first, (*it).second->x(), (*it).second->y(), (*it).second->z());
  }
  fprintf(fp, "$EndNodes\n");

  fprintf(fp, "$Elements\n");
  fprintf(fp, "%d\n", map_elem.size());
  int num_elem=1;
  for(std::map<MElement*, std::vector<unsigned long int> >::iterator it=map_elem.begin();it!=map_elem.end();++it){
    int type=whichType((*it).first);
    fprintf(fp, "%d %d %d %d %d ",num_elem,type,2,99,6); // I don't know how to access to the 3 number of tag last variable physical, geometrical region
    num_elem++;
    for(int j=0;j<(*it).second.size();j++) fprintf(fp, "%d ",(*it).second[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp, "$EndElements\n");
  fprintf(fp, "$NodeData\n");
  fprintf(fp, "1\n\"%s\"\n", name_t.c_str());
  fprintf(fp, "1\n%.16g\n", time);
  fprintf(fp, "3\n%d\n%d\n%d\n", ts, numComp, map_vertex.size());
  for(std::map<unsigned long int,SVector3>::iterator it=data.begin(); it!=data.end();++it)
    fprintf(fp, "%d %.16g %.16g %.16g\n",(*it).first,(*it).second(0),(*it).second(1),(*it).second(2));
  fprintf(fp, "$EndNodeData\n");
  fclose(fp);
  return NULL;
  }

}

/*PView *DgC0PlateSolver::buildElasticEnergyView(const std::string &postFileName)
{
  std::map<int, std::vector<double> > data;
  GaussQuadrature Integ_Bulk(GaussQuadrature::GradGrad);
  for (unsigned int i = 0; i < elasticFields.size(); ++i)
  {
    DgC0SolverField<SVector3> Field(pAssembler, LagSpace);
    IsotropicElasticTermC0Plate Eterm(Field,elasticFields[i]._E,elasticFields[i]._nu);
    DgC0BilinearTermToScalarTerm<SVector3,SVector3> Elastic_Energy_Term(Eterm);
    DgC0ScalarTermConstant One(1.0);
    for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
    {
      MElement *e=*it;
      double energ;
      double vol;
      IntPt *GP;
      int npts=Integ_Bulk.getIntPoints(e,&GP);
      Elastic_Energy_Term.get(e,npts,GP,energ);
      One.get(e,npts,GP,vol);
      std::vector<double> vec;
      vec.push_back(energ/vol);
      data[(*it)->getNum()]=vec;
    }
  }
  PView *pv = new PView (postFileName, "ElementData", pModel, data, 0.0);
  return pv;
}*/

/*PView *DgC0PlateSolver::buildVonMisesView(const std::string &postFileName)
{
  std::map<int, std::vector<double> > data;
  GaussQuadrature Integ_Bulk(GaussQuadrature::GradGrad);
  for (unsigned int i = 0; i < elasticFields.size(); ++i)
  {
    DgC0SolverField<SVector3> Field(pAssembler, LagSpace);
    IsotropicElasticTermC0Plate Eterm(Field,elasticFields[i]._E,elasticFields[i]._nu);
    DgC0BilinearTermToScalarTerm<SVector3,SVector3> Elastic_Energy_Term(Eterm);
    for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
    {
      MElement *e=*it;
      double energ;
      double vol;
      IntPt *GP;
      int npts=Integ_Bulk.getIntPoints(e,&GP);
      Elastic_Energy_Term.get(e,npts,GP,energ);
      std::vector<double> vec;
      vec.push_back(energ);
      data[(*it)->getNum()]=vec;
    }
  }
  PView *pv = new PView (postFileName, "ElementData", pModel, data, 0.0);
  return pv;
}*/

// used PView ??
void DgC0PlateSolver::buildVonMisesView(const std::string &postFileName){
  bool binary=false;               // arg of function
  std::string fileName = "stress.msh"; // arg of function
  std::string name_t = "VonMises"; //arg of function
  double time = 0.; // arg of function
  int ts =0; //arg of function What is ts;

  // used for fullDg case
  std::map<unsigned long int, MVertex*> map_vertex;
  std::set<MElement*> v; //set with element has each element has got its own dofs
  std::map<MElement*,std::vector<unsigned long int> > map_elem;

  const int numComp=1; // not arg because 1 component for VM (1 value by element)
  FILE *fp = fopen(fileName.c_str(), "w");
  // Top of file
  fprintf(fp, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
  fprintf(fp, "$Nodes\n");
  // Write node value
  // 1) Cg/Dg cases
  if(!elasticFields[0].getFormulation()){
    fprintf(fp, "%d\n", pModel->getNumVertices());
    for (unsigned int i = 0; i < elasticFields.size(); ++i)
      for (groupOfElements::vertexContainer::const_iterator it = elasticFields[i].g->vbegin(); it != elasticFields[i].g->vend(); ++it)
        fprintf(fp, "%d %.16g %.16g %.16g\n", (*it)->getNum(), (*it)->x(), (*it)->y(), (*it)->z());
  }
  else{ // 2) Full Dg Cases
    // build of set with Nodes and elements
    for (unsigned int i = 0; i < elasticFields.size(); ++i)
      for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
        v.insert(*it);
    unsigned long int num=1;
    std::vector<unsigned long int> num_vertex;
    for( std::set<MElement*>::iterator it = v.begin(); it != v.end() ;++it){
      MElement *e=*it;
      for(int j=0;j<e->getNumVertices();j++){
        map_vertex[num]=e->getVertex(j);
        num_vertex.push_back(num);
        num++;
      }
      //map_elem.insert(std::pair<int,std::vector<unsigned long int> >(e->getNum(),num_vertex));
      map_elem[e] = num_vertex;
      num_vertex.clear();
    }
    fprintf(fp, "%d\n", (int)map_vertex.size());
    for(std::map<unsigned long int, MVertex*>::iterator it = map_vertex.begin();it!=map_vertex.end();++it)
      fprintf(fp, "%d %.16g %.16g %.16g\n", (*it).first, (*it).second->x(), (*it).second->y(), (*it).second->z());
  }
  fprintf(fp, "$EndNodes\n");
  // Element data
  fprintf(fp, "$Elements\n");
  if(!elasticFields[0].getFormulation()){
    fprintf(fp, "%d\n", pModel->getNumMeshElements());
    for(int i=0;i<elasticFields.size();i++){
      for(groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it){
        int type=whichType(*it);
        fprintf(fp, "%d %d %d %d %d ",(*it)->getNum(),type,2,99,6); // I don't know how to access to the 3 number of tag last variable physical, geometrical region
      }
    }
  }
  else{
    fprintf(fp, "%d\n", map_elem.size());
    int num_elem=1;
    for(std::map<MElement*, std::vector<unsigned long int> >::iterator it=map_elem.begin();it!=map_elem.end();++it){
      int type=whichType((*it).first);
      fprintf(fp, "%d %d %d %d %d ",num_elem,type,2,99,6); // I don't know how to access to the 3 number of tag last variable physical, geometrical region
      num_elem++;
      for(int j=0;j<(*it).second.size();j++) fprintf(fp, "%d ",(*it).second[j]);
      fprintf(fp,"\n");
    }
  }
  fprintf(fp, "$EndElements\n");
  fprintf(fp, "$ElementData\n");
  fprintf(fp, "1\n\"%s\"\n", name_t.c_str());
  fprintf(fp, "1\n%.16g\n", time);
  if(!elasticFields[0].getFormulation()) fprintf(fp, "3\n%d\n%d\n%d\n", ts, numComp, pModel->getNumMeshElements());
  else fprintf(fp, "3\n%d\n%d\n%d\n", ts, numComp, map_elem.size());

  // VM value (change this ??
  GaussQuadrature Integ_Boundary(GaussQuadrature::Val);
  GaussQuadrature Integ_Bulk(GaussQuadrature::GradGrad);
  AllIPState AIPS(pModel, elasticFields, Integ_Bulk, Integ_Boundary);
  IPField<DGelasticField,DgC0FunctionSpace<SVector3> > ipf(&elasticFields,pAssembler,LagSpace,
                                                           &Integ_Bulk, &Integ_Boundary, &AIPS);
  ipf.compute1state(IPState::current);
  // Loop on element
  if(!elasticFields[0].getFormulation()){
    for (unsigned int i = 0; i < elasticFields.size(); ++i)
      for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
        fprintf(fp, "%d %.16g\n",(*it)->getNum(),ipf.getVonMises((*it),IPState::current,IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::mean));
        //fprintf(fp, "%d %.16g\n",(*it)->getNum(),10.);
  }
  else{
    int numE=0;
    for (unsigned int i = 0; i < elasticFields.size(); ++i)
      for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
        fprintf(fp, "%d %.16g\n",++numE,ipf.getVonMises((*it),IPState::current,IPField<DGelasticField,DgC0FunctionSpace<SVector3> >::mean));
        //fprintf(fp, "%d %.16g\n",++numE,10.);
  }
  fprintf(fp, "$EndElementData\n");
  fclose(fp);
};

#else
PView* DgC0PlateSolver::buildDisplacementView  (const std::string &postFileName)
{
  Msg::Error("Post-pro module not available");
  return 0;
}

PView* DgC0PlateSolver::buildElasticEnergyView(const std::string &postFileName)
{
  Msg::Error("Post-pro module not available");
  return 0;
}

// Lua interaction
#include "Bindings.h"
void DgC0PlateSolver::registerBindings(binding *b)
{
  classBinding *cb = b->addClass<DgC0PlateSolver>("DgC0PlateSolver");
  cb->setDescription("Solver class to solve plate problem with Cg/Dg or full Dg formulation");
  methodBinding *cm;
  cm = cb->addMethod("Input", &DgC0PlateSolver::readInputFile);
  cm->setArgNames("file_name",NULL);
  cm->setDescription("Use this to give a txt file with data");
  cm = cb->addMethod("CreateInterfaceElement", &DgC0PlateSolver::createInterfaceElement);
  cm->setDescription("Create all interface elements");
  cm = cb->addMethod("solve" &DgC0PlateSolver::solve);
  cm->setDescription("This function solve the linear elastic problem");
  cm = cb->addMethod("buildDisplacementView", &DgC0PlateSolver::buildDisplacementView);
  cm->setArgNames("value",NULL);
  cm->setDescription("Make a displacement view for gmsh");
  cm = cb->setConstructor<DgC0PlateSolver,int>();
  cm->setArgNames("number",NULL);
  cm->setDescription("Class of C0DgPlateSolver to solve a linear elastic plate problem with Cg/Dg or full Dg formulation");
}

#endif
