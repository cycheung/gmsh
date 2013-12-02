//////////////////////////////////////////////////////
// Templates Implementations for SystemTyped:       //
// Inclusion compilation model                      //
//                                                  //
// Damn you gcc: we want 'export' !                 //
//////////////////////////////////////////////////////

template<typename scalar>
SystemTyped<scalar>::
~SystemTyped(void){
}

template<typename scalar>
void SystemTyped<scalar>::
assemble(SolverMatrix<scalar>& A,
         SolverVector<scalar>& b,
         size_t elementId,
         const GroupOfDof& group,
         formulationPtr& term){

  const std::vector<Dof>& dof = group.getDof();
  const size_t N = group.size();

  size_t dofI;
  size_t dofJ;

  for(size_t i = 0; i < N; i++){
    dofI = dofM->getGlobalId(dof[i]);

    // If not a fixed Dof line: assemble
    if(dofI != DofManager::isFixedId()){
      for(size_t j = 0; j < N; j++){
        dofJ = dofM->getGlobalId(dof[j]);

        // If not a fixed Dof
        if(dofJ != DofManager::isFixedId())
          A.add(dofI, dofJ, (formulation->*term)(i, j, elementId));

        // If fixed Dof (for column 'dofJ'):
        //    add to right hand side (with a minus sign) !
        else
          b.add(dofI,
                -1 * dofM->getValue(dof[j]) *
                    (formulation->*term)(i, j, elementId));
      }

      b.add(dofI, formulation->rhs(i, elementId));
    }
  }
}
