///////////////////////////////////////////////////
// Templates Implementations for SystemAbstract: //
// Inclusion compilation model                   //
//                                               //
// Damn you gcc: we want 'export' !              //
///////////////////////////////////////////////////

#include <fstream>
#include "Exception.h"

template<typename scalar>
SystemAbstract<scalar>::~SystemAbstract(void){
}

template<typename scalar>
bool SystemAbstract<scalar>::isAssembled(void) const{
  return assembled;
}

template<typename scalar>
bool SystemAbstract<scalar>::isSolved(void) const{
  return solved;
}

template<typename scalar>
size_t SystemAbstract<scalar>::getSize(void) const{
  return dofM->getUnfixedDofNumber();
}

template<typename scalar>
const DofManager<scalar>& SystemAbstract<scalar>::getDofManager(void) const{
  return *dofM;
}

template<typename scalar>
const FunctionSpace& SystemAbstract<scalar>::getFunctionSpace(void) const{
  return *fs;
}

template<typename scalar>
void SystemAbstract<scalar>::constraint(const std::map<Dof, scalar>& constr){
  typename std::map<Dof, scalar>::const_iterator it  = constr.begin();
  typename std::map<Dof, scalar>::const_iterator end = constr.end();

  for(; it != end; it++)
    dofM->fixValue(it->first, it->second);
}

template<typename scalar>
void SystemAbstract<scalar>::writeMatrix(std::string fileName,
                                         std::string matrixName) const{
  std::ofstream stream;
  stream.open(fileName.c_str());
  stream << "writeMatrix not implemented for this system" << std::endl;
  stream.close();
}

template<typename scalar>
void SystemAbstract<scalar>::
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
    if(dofI != DofManager<scalar>::isFixedId()){
      for(size_t j = 0; j < N; j++){
        dofJ = dofM->getGlobalId(dof[j]);

        // If not a fixed Dof
        if(dofJ != DofManager<scalar>::isFixedId())
          A.add(dofI, dofJ, (formulation->*term)(i, j, elementId));

        // If fixed Dof (for column 'dofJ'):
        //    add to right hand side (with a minus sign) !
        else
          b.add(dofI,
                minusSign * dofM->getValue(dof[j]) *
                           (formulation->*term)(i, j, elementId));

      }

      b.add(dofI, formulation->rhs(i, elementId));
    }
  }
}
