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
          size_t order,
          scalar (*f)(fullVector<double>& xyz)){

  Basis* basis = BasisGenerator::generate(goe.get(0).getType(),
                                          0, order, "hierarchical");

  FunctionSpaceScalar                 formFS(goe, *basis);
  FormulationProjectionScalar<scalar> form(f, formFS);

  std::map<Dof, scalar> constr;
  System<scalar> projection(form);
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
          size_t order,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  Basis* basis = BasisGenerator::generate(goe.get(0).getType(),
                                          1, order, "hierarchical");

  FunctionSpaceVector                 formFS(goe, *basis);
  FormulationProjectionVector<scalar> form(f, formFS);

  std::map<Dof, scalar> constr;
  System<scalar> projection(form);
  projection.assemble();
  projection.solve();
  projection.getSolution(constr);

  sys.constraint(constr);

  delete basis;
}
