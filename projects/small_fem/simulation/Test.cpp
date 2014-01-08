#include <complex>
#include <iostream>

#include "SmallFem.h"

#include "Timer.h"

#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

#include "BasisGenerator.h"
#include "TriLagrangeBasis.h"
#include "LineNodeBasis.h"
#include "LineEdgeBasis.h"
#include "LineNedelecBasis.h"
#include "TriNodeBasis.h"
#include "QuadNedelecBasis.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationNeumann.h"
#include "FormulationSteadyWaveScalar.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "SolverMUMPS.h"

using namespace std;

complex<double> fSource(fullVector<double>& xyz){
  return complex<double>(1, 0);
}

complex<double> fWall(fullVector<double>& xyz){
  return complex<double>(0, 0);
}

int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("-msh,-o,-k");
  SmallFem::Initialize(argc, argv);

  // MPI //
  int numProcs;
  int myId;
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(numProcs != 2)
    ;//throw Exception("I just do two MPI Processes");

  // Options //
  const Options& option = SmallFem::getOptions();

  // Get Parameters //
  const double puls  = atof(option.getValue("-k")[0].c_str());
  const size_t order = atoi(option.getValue("-o")[0].c_str());

  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement* domain;
  GroupOfElement* border;
  GroupOfElement* ddm;

  if(myId == 0){
    domain = new GroupOfElement(msh.getFromPhysical(7));
    border = new GroupOfElement(msh.getFromPhysical(5));
    ddm    = new GroupOfElement(msh.getFromPhysical(4));
  }

  else{
    domain = new GroupOfElement(msh.getFromPhysical(8));
    border = new GroupOfElement(msh.getFromPhysical(6));
    ddm    = new GroupOfElement(msh.getFromPhysical(4));
  }

  // DDM Loop //
  const size_t maxIteration = 10;

  Formulation<complex<double> >* wave;
  Formulation<complex<double> >* neumann;
  System<complex<double> >*      system;

  map<Dof, complex<double> >* ddmDof;

  for(size_t k = 0; k < maxIteration; k++){
    // Formulations //
    wave =
      new FormulationSteadyWaveScalar<complex<double> > (*domain, puls, order);

    // System //
    // Init
    system = new System<complex<double> >(*wave);

    // DDM Dofs
    if(k == 0){
      set<Dof> dSet;
      set<Dof>::iterator it;
      set<Dof>::iterator end;

      system->getFunctionSpace().getKeys(*ddm, dSet);

      it     = dSet.begin();
      end    = dSet.end();
      ddmDof = new std::map<Dof, complex<double> >;

      for(; it != end; it++)
        ddmDof->insert(pair<Dof, complex<double> >(*it, 0));
    }

    if(myId == 0){// && (k == 0 || k == maxIteration - 1)){
      map<Dof, complex<double> >::iterator it  = ddmDof->begin();
      map<Dof, complex<double> >::iterator end = ddmDof->end();

      for(; it != end; it++)
        cout << it->first.toString() << ": " << it->second << endl;
      cout << " --- " << endl;
    }

    // Constraint
    if(myId == 0)
      SystemHelper<complex<double> >::dirichlet(*system, *border, fSource);
    else
      SystemHelper<complex<double> >::dirichlet(*system, *border, fWall);

    // Assemble
    system->assemble();

    // Assemble Neumann term
    neumann = new FormulationNeumann(*ddm, puls, order);
    system->addBorderTerm(*neumann);

    // Solve
    system->solve();

    // Get DDM Solution //
    system->getSolution(*ddmDof, 0);

    // Write Solution //
    stringstream feSolName;
    FEMSolution<complex<double> > feSol;
    system->getSolution(feSol);

    if(myId == 0)
      feSolName << "proc0_it" << k;
    else
      feSolName << "proc1_it" << k;

    feSol.write(feSolName.str());

    // Clean //
    delete neumann;
    delete wave;
    delete system;
  }

  // Finalize //
  delete domain;
  delete border;
  delete ddm;
  delete ddmDof;

  SmallFem::Finalize();
}
