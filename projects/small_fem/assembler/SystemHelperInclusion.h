/////////////////////////////////////////////////
// Templates Implementations for SystemHelper: //
// Inclusion compilation model                 //
//                                             //
// Damn you gcc: we want 'export' !            //
/////////////////////////////////////////////////

#include "System.h"
#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"


template<typename scalar>
SystemHelper<scalar>::SystemHelper(void){
}

template<typename scalar>
SystemHelper<scalar>::~SystemHelper(void){
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          GroupOfElement& goe,
          scalar (*f)(fullVector<double>& xyz)){

  const FunctionSpace& fs = sys.getFunctionSpace();

  Basis* basis = BasisGenerator::generate(goe.get(0).getType(),
                                          fs.getBasis(0).getType(),
                                          fs.getBasis(0).getOrder(),
                                          "hierarchical");

  FunctionSpaceScalar         formFS(goe, *basis);
  FormulationProjectionScalar form(f, formFS);

  std::map<Dof, double> constr;
  System projection(form);
  projection.assemble();
  projection.solve();
  projection.getSolution(constr);

  sys.constraint(constr);

  delete basis;
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          GroupOfElement& goe,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  const FunctionSpace& fs = sys.getFunctionSpace();

  Basis* basis = BasisGenerator::generate(goe.get(0).getType(),
                                          fs.getBasis(0).getType(),
                                          fs.getBasis(0).getOrder(),
                                          "hierarchical");

  FunctionSpaceVector         formFS(goe, *basis);
  FormulationProjectionVector form(f, formFS);

  std::map<Dof, double> constr;
  System projection(form);
  projection.assemble();
  projection.solve();
  projection.getSolution(constr);

  sys.constraint(constr);

  delete basis;
}
