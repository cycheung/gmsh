#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "Mesh.h"
#include "System.h"
#include "fullMatrix.h"
#include "FEMSolution.h"
#include "Interpolator.h"
#include "BasisGenerator.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"

#include "SmallFem.h"

using namespace std;

//////////////////////////////////////////////
// Functions to project (scalar and vector) //
//////////////////////////////////////////////

double fScal(fullVector<double>& xyz){
  return
    sin(10 * xyz(0)) +
    sin(10 * xyz(1)) +
    sin(10 * xyz(2));
}

fullVector<double> fVect(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = sin(10 * xyz(0));
  res(1) = sin(10 * xyz(1));
  res(2) = sin(10 * xyz(2));

  return res;
}

/////////////////////////////////////
// Prototypes of helping functions //
/////////////////////////////////////

// Analytical solution (scalar and vector) //
void ana(double (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<double>& eval);

void ana(fullVector<double> (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<double>& eval);

// FEM solution (scalar and vector) //
void fem(double (*f)(fullVector<double>& xyz),
         GroupOfElement& domain,
         size_t order,
         const fullMatrix<double>& point,
         fullMatrix<double>& sol,
         bool nopos);

void fem(fullVector<double> (*f)(fullVector<double>& xyz),
         GroupOfElement& domain,
         size_t order,
         const fullMatrix<double>& point,
         fullMatrix<double>& sol,
         bool nopos);

// L2 norms and error (|A|, |A - B| and |A - B| / |A|) //
double l2Norm(const fullMatrix<double>& val);
double l2Norm(const fullMatrix<double>& valA, const fullMatrix<double>& valB);
double l2Error(const fullMatrix<double>& fem, const fullMatrix<double>& ana);

// Write Octave file //
void write(bool isScalar, const fullMatrix<double>& l2, string name);

////////////////////////////////////
// SmallFem Computations and main //
////////////////////////////////////

void compute(const Options& option){
  // Get Visu Mesh and its vertex coordinates //
  cout << "## Reference Mesh" << endl << flush;
  Mesh           visuMsh(option.getValue("-ref")[0]);
  GroupOfElement visu = visuMsh.getFromPhysical(7);

  fullMatrix<double> point;
  visuMsh.getAllVertexCoordinate(point);

  // Get FEM Orders //
  const size_t nOrder = option.getValue("-o").size();
  vector<int>   order(nOrder);

  for(size_t i = 0; i < nOrder; i++)
    order[i] = atoi(option.getValue("-o")[i].c_str());

  // Get FEM Meshes //
  const size_t  nMesh = option.getValue("-msh").size();
  vector<string> mesh(nMesh);

  for(size_t i = 0; i < nMesh; i++)
    mesh[i] = option.getValue("-msh")[i];

  // Post Processing ? //
  bool nopos = (option.getValue("-nopos").size() != 0);

  // Scalar or Vector //
  bool isScalar = (option.getValue("-type")[0].compare("scalar") == 0);

  // Real Solutions //
  cout << "## Real Solution" << endl << flush;
  fullMatrix<double> realSol;

  if(isScalar)
    ana(fScal, point, realSol);
  else
    ana(fVect, point, realSol);

  // FEM Solution & L2 Error //
  cout << "## FEM Solutions" << endl << flush;
  fullMatrix<double> femSol;
  fullMatrix<double> l2(nOrder, nMesh);

  // Iterate on Meshes
  for(size_t i = 0; i < nMesh; i++){
    cout << " ** Mesh: " << mesh[i] << endl << flush;
    Mesh           msh(mesh[i]);
    GroupOfElement domain = msh.getFromPhysical(7);

    // Iterate on Orders
    for(size_t j = 0; j < nOrder; j++){
      cout << "  -- Order " << order[j] << ": " << flush;

      if(isScalar)
        fem(fScal, domain, order[j], point, femSol, nopos);
      else
        fem(fVect, domain, order[j], point, femSol, nopos);

      l2(j, i) = l2Error(femSol, realSol);
      cout << l2(j, i) << endl;
    }
  }

  // Display //
  cout << "## L2 Error" << endl << flush;
  write(isScalar, l2, option.getValue("-name")[0]);
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-type,-ref,-nopos,-name");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}

///////////////////////////////////////
// Implentation of helping functions //
///////////////////////////////////////

void ana(double (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<double>& eval){
  // Alloc eval for Scalar Values //
  const size_t nPoint = point.size1();
  eval.resize(nPoint, 1);

  // Loop on point and evaluate f //
  fullVector<double> xyz(3);
  for(size_t i = 0; i < nPoint; i++){
    xyz(0) = point(i, 0);
    xyz(1) = point(i, 1);
    xyz(2) = point(i, 2);

    eval(i, 0) = f(xyz);
  }
}

void ana(fullVector<double> (*f)(fullVector<double>& xyz),
         const fullMatrix<double>& point,
         fullMatrix<double>& eval){
  // Alloc eval for Vectorial Values //
  const size_t nPoint = point.size1();
  eval.resize(nPoint, 3);

  // Loop on point and evaluate f //
  fullVector<double> xyz(3);
  fullVector<double> tmp(3);
  for(size_t i = 0; i < nPoint; i++){
    xyz(0) = point(i, 0);
    xyz(1) = point(i, 1);
    xyz(2) = point(i, 2);

    tmp = f(xyz);

    eval(i, 0) = tmp(0);
    eval(i, 1) = tmp(1);
    eval(i, 2) = tmp(2);
  }
}

void fem(double (*f)(fullVector<double>& xyz),
         GroupOfElement& domain,
         size_t order,
         const fullMatrix<double>& point,
         fullMatrix<double>& sol,
         bool nopos){
  // Basis//
  stringstream stream;
  Basis* basis  = BasisGenerator::generate(domain.get(0).getType(),
                                           0, order, "hierarchical");
  // Projection //
  FunctionSpaceScalar fSpace(domain, *basis);
  FormulationProjectionScalar projection(f, fSpace);
  System sysProj(projection);

  // Assemble and Solve //
  sysProj.assemble();
  sysProj.solve();

  // Interpolate on Ref Points //
  fullVector<double> sysSol;
  sysProj.getSolution(sysSol);

  Interpolator::interpolate(sysProj.getFunctionSpace(),
                            sysProj.getDofManager(),
                            sysSol,
                            point,
                            sol);
  // Post-processing //
  if(!nopos){
    FEMSolution feSol;
    stream << "projection_Mesh" << domain.getNumber() << "_Order" << order;

    sysProj.getSolution(feSol);
    feSol.write(stream.str());
  }
}

void fem(fullVector<double> (*f)(fullVector<double>& xyz),
         GroupOfElement& domain,
         size_t order,
         const fullMatrix<double>& point,
         fullMatrix<double>& sol,
         bool nopos){
  // Basis//
  stringstream stream;
  Basis* basis  = BasisGenerator::generate(domain.get(0).getType(),
                                           1, order, "hierarchical");
  // Projection //
  FunctionSpaceVector fSpace(domain, *basis);
  FormulationProjectionVector projection(f, fSpace);
  System sysProj(projection);

  // Assemble and Solve //
  sysProj.assemble();
  sysProj.solve();

  // Interpolate on Ref Points //
  fullVector<double> sysSol;
  sysProj.getSolution(sysSol);

  Interpolator::interpolate(sysProj.getFunctionSpace(),
                            sysProj.getDofManager(),
                            sysSol,
                            point,
                            sol);
  // Post-processing //
  if(!nopos){
    FEMSolution feSol;
    stream << "projection_Mesh" << domain.getNumber() << "_Order" << order;

    sysProj.getSolution(feSol);
    feSol.write(stream.str());
  }
}

double l2Norm(const fullMatrix<double>& val){
  const size_t nPoint = val.size1();
  const size_t    dim = val.size2();

  double norm = 0;
  double modSquare = 0;

  for(size_t i = 0; i < nPoint; i++){
    modSquare = 0;

    for(size_t j = 0; j < dim; j++)
      modSquare += val(i, j) * val(i, j);

    norm += modSquare;
  }

  return norm = sqrt(norm);
}

double l2Norm(const fullMatrix<double>& valA, const fullMatrix<double>& valB){
  const size_t nPoint = valA.size1();
  const size_t    dim = valA.size2();

  double norm = 0;
  double modSquare = 0;

  for(size_t i = 0; i < nPoint; i++){
    modSquare = 0;

    for(size_t j = 0; j < dim; j++)
      modSquare += (valA(i, j) - valB(i, j)) * (valA(i, j) - valB(i, j));

    norm += modSquare;
  }

  return norm = sqrt(norm);
}

double l2Error(const fullMatrix<double>& fem, const fullMatrix<double>& ana){
  double anaNorm = l2Norm(ana);
  double femNorm = l2Norm(fem, ana);

  return femNorm / anaNorm;
}

void write(bool isScalar, const fullMatrix<double>& l2, string name){
  // Stream
  ofstream stream;
  string   fileName;
  if(isScalar)
    fileName = name + "Node.m";
  else
    fileName = name + "Edge.m";

  stream.open(fileName.c_str());

  // Matrix data
  const size_t l2Row      = l2.size1();
  const size_t l2ColMinus = l2.size2() - 1;

  // Clean Octave
  stream << "close all;" << endl
         << "clear all;" << endl
         << endl;

  // Mesh (Assuming uniform refinement)
  stream << "h = [1, ";
  for(size_t i = 1; i < l2ColMinus; i++)
    stream << 1 / pow(2, i) << ", ";

  stream << 1 / pow(2, l2ColMinus) << "];" << endl;

  // Order (Assuming uniform refinement)
  stream << "p = [1:" << l2Row << "];" << endl
         << endl;

  // Matrix of l2 error (l2[Order][Mesh])
  stream << "l2 = ..." << endl
         << "    [..." << endl;

  for(size_t i = 0; i < l2Row; i++){
    stream << "        ";

    for(size_t j = 0; j < l2ColMinus; j++)
      stream << scientific << showpos
             << l2(i, j) << " , ";

    stream << scientific << showpos
           << l2(i, l2ColMinus) << " ; ..." << endl;
  }

  stream << "    ];" << endl
         << endl;

  // Delta
  stream << "P = size(p, 2);"                                   << endl
         << "H = size(h, 2);"                                   << endl
         << endl
         << "delta = zeros(P, H - 1);"                          << endl
         << endl
         << "for i = 1:H-1"                                     << endl
         << "  delta(:, i) = ..."                               << endl
         << "    (log10(l2(:, i + 1)) - log10(l2(:, i))) / ..." << endl
         << "    (log10(1/h(i + 1))   - log10(1/h(i)));"        << endl
         << "end"                                               << endl
         << endl
         << "delta"                                             << endl
         << endl;

  // Plot
  stream << "figure;"                  << endl
         << "loglog(1./h, l2', '-*');" << endl
         << "grid;"
         << endl;

  // Title
  stream << "title(" << "'" << name << ": ";
  if(isScalar)
    stream << "Nodal";
  else
    stream << "Edge";

  stream << "');" << endl
         << endl;

  // Axis
  stream << "xlabel('1/h [-]');"      << endl
         << "ylabel('L2 Error [-]');" << endl;
  // Close
  stream.close();
}
