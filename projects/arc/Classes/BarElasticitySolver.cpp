//
// Description : BAR elasticity solver, element space function enriched on tagged vertex
//
//
// Author:  <Eric Bechet>::<Boris Sedji>,  01/2010
//
// Copyright: See COPYING file that comes with this distribution
//
//



#include <string.h>
#include "GmshConfig.h"
#include "elasticitySolver.h"
#include "linearSystemCSR.h"
#include "linearSystemPETSc.h"
#include "linearSystemGMM.h"
#include "Numeric.h"
#include "GModel.h"
#include "functionSpace.h"
#include "terms.h"
#include "solverAlgorithms.h"
#include "quadratureRules.h"
#include "solverField.h"
#if defined(HAVE_POST)
#include "PView.h"
#include "PViewData.h"
#include "MPoint.h"
#endif


#include "OctreeSolver.h"
#include "DILevelset.h"
#include "NodeEnrichedFS.cpp"
#include "gLevelSetMesh.h"

#include "BarElasticitySolver.h"

#include "groupOfElements.h"
#include "FuncGradDisc.h"
#include "FuncHeaviside.h"

#include "xFemFunctionSpace.cpp"

#include "BeamParam.h"
#include "BeamTerm.h"
#include "BeamLin.h"

BarElasticitySolver::~BarElasticitySolver()
{
  delete _funcEnrichment;
}


void BarElasticitySolver::setMesh(const std::string &meshFileName)
{
  pModel = new GModel();
  pModel->readMSH(meshFileName.c_str());
   _dim = pModel->getNumRegions() ? 3 : 2;
}

void BarElasticitySolver::solve(){

#if defined(HAVE_TAUCS)
  linearSystemCSRTaucs<double> *lsys = new linearSystemCSRTaucs<double>;
#elif defined(HAVE_PETSC)
  linearSystemPETSc<double> *lsys = new linearSystemPETSc<double>;
#else
  linearSystemGmm<double> *lsys = new linearSystemGmm<double>;
  lsys->setNoisy(2);
#endif
  if (pAssembler) delete pAssembler;
  pAssembler = new dofManager<double>(lsys);

  // determine all enriched nodes and save in map member _TagEnrichedVertex

//  for (int i = 0; i < elasticFields.size(); ++i){
//    groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin();
//    for ( ; it != elasticFields[i].g->end(); ++it)
//    {
//      MElement *e = *it;
//      if (e->getParent()) // if element got parents
//      {
//        for (int k = 0; k < e->getParent()->getNumVertices(); ++k)
//        {  // for all vertices in the element parent
//          _TagEnrichedVertex.insert(e->getParent()->getVertex(k)->getNum());
//        }
//      }
//    }
//  }

  _TagEnrichedVertex.clear();


  // enriched composant
  _EnrichComp.insert(0);
  _EnrichComp.insert(1);
  //_EnrichComp.push_back(2);

  // level set definition (in .dat ??)
  double a(0.), b(1.), c(0.), d(-4.7);
  int n(1); // pour level set
  gLevelset *ls = new gLevelsetPlane(a, b, c, d, n);
  _funcEnrichment = new FuncHeaviside(ls);


  // space function definition
  FunctionSpace<SVector3> *NLagSpace;
  if (LagSpace) delete LagSpace;
  if (_dim==3) NLagSpace=new VectorLagrangeFunctionSpace(_tag);
  if (_dim==2) NLagSpace=new VectorLagrangeFunctionSpace(_tag,VectorLagrangeFunctionSpace::VECTOR_X,VectorLagrangeFunctionSpace::VECTOR_Y);


  // space function definition //
  // Lagrange space with normal dofs
  FunctionSpace<SVector3> *LagrangeSpace;
  if (_dim==3) LagrangeSpace=new VectorLagrangeFunctionSpace(_tag);
  if (_dim==2) LagrangeSpace=new VectorLagrangeFunctionSpace(_tag,VectorLagrangeFunctionSpace::VECTOR_X,VectorLagrangeFunctionSpace::VECTOR_Y);

  // Filtered space
  FilterNodeEnriched *filter = new FilterNodeEnriched(&_TagEnrichedVertex , &_EnrichComp);
  FilteredFS<SVector3> *LagrangeFiltered = new FilteredFS<SVector3>(LagrangeSpace,filter);

  // xfem filtered space
  xFemFS<SVector3> *xFemLagrange = new xFemFS<SVector3>(LagrangeFiltered,_funcEnrichment);

  // Composite space : xfem + Lagrange
  if (LagSpace) delete LagSpace;

  LagSpace = new CompositeFunctionSpace<SVector3>(*LagrangeSpace, *xFemLagrange);


  // Beams Definition
  std::vector <SPoint2> Points;
  SPoint2 P(5,90);
  SPoint2 Q(5,10);
  Points.push_back(P);
  Points.push_back(Q);
  BeamLin Beam1(Points);

  for (int i = 0 ; i < 1000 ;i++)  _Beams.push_back(Beam1);

  std::cout <<  "Dirichlet BC"<< std::endl;
  for (unsigned int i = 0; i < allDirichlet.size(); i++)
  {
    FilterDofComponent filter(allDirichlet[i]._comp);
    FixNodalDofs(*LagSpace,allDirichlet[i].g->begin(),allDirichlet[i].g->end(),*pAssembler,allDirichlet[i]._f,filter);
  }

  // we number the dofs : when a dof is numbered, it cannot be numbered
  // again with another number.
  for (unsigned int i = 0; i < elasticFields.size(); ++i)
  {
    NumberDofs(*LagSpace, elasticFields[i].g->begin(), elasticFields[i].g->end(),*pAssembler);
  }

  // Now we start the assembly process
  // First build the force vector

  GaussQuadrature Integ_Boundary(GaussQuadrature::Val);
  std::cout <<  "Neumann BC"<< std::endl;
  for (unsigned int i = 0; i < allNeumann.size(); i++)
  {
    LoadTerm<SVector3> Lterm(*LagSpace,allNeumann[i]._f);
    Assemble(Lterm,*LagSpace,allNeumann[i].g->begin(),allNeumann[i].g->end(),Integ_Boundary,*pAssembler);
  }

// bulk material law

  GaussQuadrature Integ_Bulk(GaussQuadrature::GradGrad);
  for (unsigned int i = 0; i < elasticFields.size(); i++)
  {
    IsotropicElasticTerm Eterm(*LagSpace,elasticFields[i]._E,elasticFields[i]._nu);
//    LaplaceTerm<SVector3,SVector3> Eterm(*LagSpace);
    Assemble(Eterm,*LagSpace,elasticFields[i].g->begin(),elasticFields[i].g->end(),Integ_Bulk,*pAssembler);
  }

  // Beams Integration
  BeamsIntegration();


  printf("-- done assembling!\n");
  std::cout<<"Dof Number : " << pAssembler->sizeOfR() <<"\n";
  lsys->systemSolve();
  printf("-- done solving!\n");
  double energ=0;
  for (unsigned int i = 0; i < elasticFields.size(); i++)
  {
    SolverField<SVector3> Field(pAssembler, LagSpace);
    IsotropicElasticTerm Eterm(Field,elasticFields[i]._E,elasticFields[i]._nu);
    BilinearTermToScalarTerm<SVector3,SVector3> Elastic_Energy_Term(Eterm);
    Assemble(Elastic_Energy_Term,elasticFields[i].g->begin(),elasticFields[i].g->end(),Integ_Bulk,energ);
  }
  printf("elastic energy=%f\n",energ);
}


void BarElasticitySolver::BeamsIntegration()
{
  GaussQuadrature Integ_Bulk(GaussQuadrature::GradGrad);
  for (unsigned int i = 0; i < elasticFields.size(); ++i)
  {
    for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
    {
      MElement *e=*it;
      if (e->getParent()) e=e->getParent();
      for (int j = 0; j < _Beams.size();j++)
      {
        BeamParam * BeamElem = _Beams[j].cut(e) ;
        if (BeamElem != NULL)
        {
        BeamTerm BTerm(LagSpace,e,BeamElem);
        fullMatrix<double> localMatrix;
        std::vector<Dof> R;
        R.clear();
        //IntPt *GP;
        //int npts=Integ_Bulk.getIntPoints(e,&GP);
        //BTerm.get(npts,GP,localMatrix);
        BTerm.get(localMatrix);
        LagSpace->getKeys(e,R);
        pAssembler->assemble(R, localMatrix);
        }
      }
    }
  }
}

PView* BarElasticitySolver::buildDisplacementView (const std::string &postFileName)
{

  std::set<MVertex*> v;
  for (unsigned int i = 0; i < elasticFields.size(); ++i)
  {
    for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it)
    {
      MElement *e=*it;
      if (e->getParent()) e=e->getParent();
      for (int j = 0; j < e->getNumVertices(); ++j) v.insert(e->getVertex(j));
    }
  }
  std::map<int, std::vector<double> > data;
  SolverField<SVector3> Field(pAssembler, LagSpace);
 for ( std::set<MVertex*>::iterator it = v.begin(); it != v.end(); ++it)
  {
    SVector3 val;
    MPoint p(*it);
    Field.f(&p,0,0,0,val);
    std::vector<double> vec(3);vec[0]=val(0);vec[1]=val(1);vec[2]=val(2);
    //std::cout<<"ver num : " << (*it)->getNum() <<"\n" ;
    //std::cout<<"( "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" )\n";
    data[(*it)->getNum()]=vec;
  }
  PView *pv = new PView (postFileName, "NodeData", pModel, data, 0.0);
  return pv;

}
