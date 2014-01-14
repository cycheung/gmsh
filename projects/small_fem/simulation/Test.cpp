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
#include "FormulationEMDA.h"
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
    throw Exception("I just do two MPI Processes");

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
  Formulation<complex<double> >* emda;
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

    // Constraint
    if(myId == 0)
      SystemHelper<complex<double> >::dirichlet(*system, *border, fSource);
    else
      SystemHelper<complex<double> >::dirichlet(*system, *border, fWall);

    // Assemble
    system->assemble();
    /*
    if(myId == 1){// && (k == 0 || k == maxIteration - 1)){
      map<Dof, complex<double> >::iterator it  = ddmDof->begin();
      map<Dof, complex<double> >::iterator end = ddmDof->end();

      for(; it != end; it++)
        cout << "g" << myId << ": " << it->first.toString()
             << ": " << it->second << endl;
      cout << " --- " << endl;
    }
    */
    // Neumann and EMDA terms
    neumann = new FormulationNeumann(*ddm, puls, order);
    emda    = new FormulationEMDA
      (static_cast<const FunctionSpaceScalar&>(neumann->fs()), *ddmDof);

    system->addBorderTerm(*neumann);
    system->addBorderTerm(*emda);

    // Solve
    system->solve();

    // Get DDM Solution //
    map<Dof, complex<double> > oldDdmDof = *ddmDof;
    system->getSolution(*ddmDof, 0);
    /*
    if(myId == 1){// && (k == 0 || k == maxIteration - 1)){
      map<Dof, complex<double> >::iterator it  = ddmDof->begin();
      map<Dof, complex<double> >::iterator end = ddmDof->end();

      for(; it != end; it++)
        cout << "u" << myId << ": " << it->first.toString()
             << ": " << it->second << endl;
      cout << " --- " << endl;
    }
    */
    // Update DDM //
    // Upade my Values
    map<Dof, complex<double> >::iterator it;
    map<Dof, complex<double> >::iterator end;

    map<Dof, complex<double> >::iterator it2;
    map<Dof, complex<double> >::iterator end2;

    it  = ddmDof->begin();
    end = ddmDof->end();

    it2  = oldDdmDof.begin();
    end2 = oldDdmDof.end();

    for(; it != end; it++, it2++){
      if(it->first == it2->first){
        // g_new = 2*j*k*u + g_old
        it->second =
          (it->second * (complex<double>(0, 2 * puls))) + it2->second;
      }

      else
        throw Exception("Snif");
    }
    // Serialize my Values
    const size_t             ddmDofSize = ddmDof->size();
    vector<int>              ddmDofEntity(ddmDofSize, 0);
    vector<int>              ddmDofType(ddmDofSize, 0);
    vector<complex<double> > ddmDofValue(ddmDofSize, 0);

    vector<int>              incomingDdmDofEntity(ddmDofSize, 0);
    vector<int>              incomingDdmDofType(ddmDofSize, 0);
    vector<complex<double> > incomingDdmDofValue(ddmDofSize, 0);

    it  = ddmDof->begin();
    end = ddmDof->end();

    for(size_t i = 0; it != end; i++, it++){
      ddmDofEntity[i] = (int)(it->first.getEntity());
      ddmDofType[i]   = (int)(it->first.getType());
      ddmDofValue[i]  = it->second;
    }

    // Exchange
    if(myId == 0){
      // Send to 1
      MPI_Status status;

      MPI_Ssend((void*)(ddmDofEntity.data()), ddmDofSize,
                MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Ssend((void*)(ddmDofType.data()), ddmDofSize,
                MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Ssend((void*)(ddmDofValue.data()), ddmDofSize,
                MPI::DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD);

      // Recv from 1
      MPI_Recv((void*)(incomingDdmDofEntity.data()), ddmDofSize,
               MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv((void*)(incomingDdmDofType.data()), ddmDofSize,
               MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv((void*)(incomingDdmDofValue.data()), ddmDofSize,
               MPI::DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD, &status);
    }

    if(myId == 1){
      // Recv from 0
      MPI_Status status;

      MPI_Recv((void*)(incomingDdmDofEntity.data()), ddmDofSize,
               MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv((void*)(incomingDdmDofType.data()), ddmDofSize,
               MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv((void*)(incomingDdmDofValue.data()), ddmDofSize,
               MPI::DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &status);

      // Send to 0
      MPI_Ssend((void*)(ddmDofEntity.data()), ddmDofSize,
                MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Ssend((void*)(ddmDofType.data()), ddmDofSize,
                MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Ssend((void*)(ddmDofValue.data()), ddmDofSize,
                MPI::DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
    }

    // Update ddmDof
    it  = ddmDof->begin();
    end = ddmDof->end();

    it2  = oldDdmDof.begin();
    end2 = oldDdmDof.end();

    for(size_t i = 0; it != end; it++, it2++, i++){
      if((int)(it->first.getType()) != incomingDdmDofType[i])
        throw Exception("Snif");

      if((int)(it->first.getEntity()) != incomingDdmDofEntity[i])
        throw Exception("Snif");

      if((it->first.getType()) != it2->first.getType())
        throw Exception("Snif");

      if((it->first.getEntity()) != it2->first.getEntity())
        throw Exception("Snif");

      it->second = it2->second + incomingDdmDofValue[i];
    }

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
    delete emda;
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
