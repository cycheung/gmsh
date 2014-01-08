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

  // Get Function Space for Projection (formFS) //
  const FunctionSpace& fs = sys.getFunctionSpace();

  Basis* basis = BasisGenerator::generate(goe.get(0).getType(),
                                          fs.getBasis(0).getType(),
                                          fs.getBasis(0).getOrder(),
                                          "hierarchical");

  FunctionSpaceScalar formFS(goe, *basis);

  // Solve Projection //
  FormulationProjectionScalar<scalar> form(f, formFS);

  System<scalar> projection(form);
  projection.assemble();
  projection.solve();

  // Map of Dofs //
  const std::set<Dof>& dof = formFS.getAllDofs();
  std::set<Dof>::iterator it  = dof.begin();
  std::set<Dof>::iterator end = dof.end();

  std::map<Dof, scalar> constr;
  for(; it != end; it++)
    constr.insert(std::pair<Dof, scalar>(*it, 0));

  // Get Solution and Dirichlet Constraint //
  projection.getSolution(constr, 0);
  sys.constraint(constr);

  delete basis;
}

template<typename scalar>
void SystemHelper<scalar>::
dirichlet(SystemAbstract<scalar>& sys,
          GroupOfElement& goe,
          fullVector<scalar> (*f)(fullVector<double>& xyz)){

  // Get Function Space for Projection (formFS) //
  const FunctionSpace& fs = sys.getFunctionSpace();

  Basis* basis = BasisGenerator::generate(goe.get(0).getType(),
                                          fs.getBasis(0).getType(),
                                          fs.getBasis(0).getOrder(),
                                          "hierarchical");

  FunctionSpaceVector formFS(goe, *basis);

  // Solve Projection //
  FormulationProjectionVector<scalar> form(f, formFS);

  System<scalar> projection(form);
  projection.assemble();
  projection.solve();

  // Map of Dofs //
  const std::set<Dof>& dof = formFS.getAllDofs();
  std::set<Dof>::iterator it  = dof.begin();
  std::set<Dof>::iterator end = dof.end();

  std::map<Dof, scalar> constr;
  for(; it != end; it++)
    constr.insert(std::pair<Dof, scalar>(*it, 0));

  // Get Solution and Dirichlet Constraint //
  projection.getSolution(constr, 0);
  sys.constraint(constr);

  delete basis;
}
