#include <cmath>
#include <complex>
#include <iostream>

#include "SmallFem.h"

#include "Timer.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationNeumann.h"
#include "FormulationSteadyWaveScalar.h"

#include "FormulationEMDA.h"
#include "FormulationUpdateEMDA.h"

#include "FormulationOO2.h"
#include "FormulationUpdateOO2.h"

using namespace std;

typedef std::complex<double> Complex;

Complex fSource(fullVector<double>& xyz){
  return Complex(1, 0);
}

void initMap(System<Complex>& system,
             GroupOfElement& goe,
             map<Dof, Complex>& map){

  set<Dof> dSet;
  system.getFunctionSpace().getKeys(goe, dSet);

  set<Dof>::iterator it  = dSet.begin();
  set<Dof>::iterator end = dSet.end();

  for(; it != end; it++)
    map.insert(pair<Dof, Complex>(*it, 0));
}

void displaySolution(map<Dof, Complex>& solution,
                     int id, size_t step, string name){
  int myId;
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(myId == id){
    cout << "After Iteration: " << step + 1 << endl;
    map<Dof, Complex>::iterator it  = solution.begin();
    map<Dof, Complex>::iterator end = solution.end();

    for(; it != end; it++)
      cout << name << myId << ": " << it->first.toString()
           << ": " << it->second << endl;
    cout << " --- " << endl;
  }
}

void serialize(const map<Dof, Complex>& map,
               vector<int>& entity, vector<int>& type, vector<Complex>& value){

  std::map<Dof, Complex>::const_iterator it  = map.begin();
  std::map<Dof, Complex>::const_iterator end = map.end();

  for(size_t i = 0; it != end; i++, it++){
    entity[i] = (int)(it->first.getEntity());
    type[i]   = (int)(it->first.getType());
    value[i]  = it->second;
  }
}

void unserialize(map<Dof, Complex>& map,
                 const vector<int>& entity, const vector<int>& type,
                 const vector<Complex>& value){

  std::map<Dof, Complex>::iterator it  = map.begin();
  std::map<Dof, Complex>::iterator end = map.end();

  for(size_t i = 0; it != end; it++, i++){
    if((int)(it->first.getType()) != type[i])
      throw Exception("Snif");

    if((int)(it->first.getEntity()) != entity[i])
      throw Exception("Snif");

    it->second = value[i];
  }
}

void exchange(int myId,
              vector<int>& outEntity, vector<int>& outType,
              vector<Complex>& outValue,
              vector<int>& inEntity, vector<int>& inType,
              vector<Complex>& inValue){

  MPI_Status s;
  size_t     size = outEntity.size();

  if(myId == 0){
    // Send to 1
    MPI_Ssend((void*)(outEntity.data()), size, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outType.data()),   size, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outValue.data()),  size,
              MPI::DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD);

    // Recv from 1
    MPI_Recv((void*)(inEntity.data()), size, MPI_INT, 1, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inType.data()),   size, MPI_INT, 1, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inValue.data()),  size,
             MPI::DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD, &s);
  }

  if(myId == 1){
    // Recv from 0
    MPI_Recv((void*)(inEntity.data()), size, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inType.data()),   size, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
    MPI_Recv((void*)(inValue.data()),  size,
             MPI::DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &s);

    // Send to 0
    MPI_Ssend((void*)(outEntity.data()), size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outType.data()),   size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Ssend((void*)(outValue.data()),  size,
              MPI::DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
  }
}

void compute(const Options& option){
  // MPI //
  int numProcs;
  int myId;
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myId);

  if(numProcs != 2)
    ;//throw Exception("I just do two MPI Processes");

  // Get Parameters //
  const string ddmType = option.getValue("-ddm")[1];
  const double k       = atof(option.getValue("-k")[1].c_str());
  const size_t order   = atoi(option.getValue("-o")[1].c_str());
  const size_t maxIt   = atoi(option.getValue("-max")[1].c_str());

  // DDM Formulations //
  const string emdaType("emda");
  const string oo2Type("oo2");

  // Variables
  const double Pi = std::atan(1.0) * 4;
  double lc       = 0;
  double chi      = 0;
  Complex ooA     = 0;
  Complex ooB     = 0;

  // EMDA Stuff
  if(ddmType == emdaType)
    chi = atof(option.getValue("-chi")[1].c_str());

  // OO2 Stuff
  if(ddmType == oo2Type){
    lc = atof(option.getValue("-lc")[1].c_str());

    double ooXsiMin = 0;
    double ooXsiMax = Pi / lc;
    double ooDeltaK = Pi / .06;

    double tmp0 =
      (k * k - ooXsiMin * ooXsiMin) * (k * k - (k - ooDeltaK) * (k - ooDeltaK));

    double tmp1 =
      (ooXsiMax * ooXsiMax - k * k) * ((k + ooDeltaK) * (k + ooDeltaK) - k * k);

    Complex ooAlpha = std::pow(Complex(tmp0, 0), 0.25) * Complex(0, 1);
    Complex ooBeta  = std::pow(Complex(tmp1, 0), 0.25);

    ooA = -(ooAlpha * ooBeta - k * k) / (ooAlpha + ooBeta);
    ooB = Complex(-1, 0) / (ooAlpha + ooBeta);
  }

  // Get Domains //
  Mesh msh(option.getValue("-msh")[1]);
  GroupOfElement* domain;
  GroupOfElement* source;
  GroupOfElement* infinity;
  GroupOfElement* ddmBorder;

  if(myId == 0){
    domain    = new GroupOfElement(msh.getFromPhysical(7));
    source    = new GroupOfElement(msh.getFromPhysical(5));
    infinity  = new GroupOfElement(msh.getFromPhysical(61));
    ddmBorder = new GroupOfElement(msh.getFromPhysical(4));
  }

  else{
    domain    = new GroupOfElement(msh.getFromPhysical(8));
    source    = NULL;
    infinity  = new GroupOfElement(msh.getFromPhysical(62));
    ddmBorder = new GroupOfElement(msh.getFromPhysical(4));
  }

  // Formulation Pointers //
  Formulation<Complex>* wave;
  Formulation<Complex>* ddm;
  Formulation<Complex>* upDdm;
  Formulation<Complex>* neumann;

  // System Pointers //
  System<Complex>* system;
  System<Complex>* update;

  // Solution Maps Pointers //
  map<Dof, Complex>* ddmG     = NULL;
  map<Dof, Complex>* solution = NULL;

  // MPI Buffers Pointers //
  vector<int>*     outEntity = NULL;
  vector<int>*     outType   = NULL;
  vector<Complex>* outValue  = NULL;

  vector<int>*     inEntity  = NULL;
  vector<int>*     inType    = NULL;
  vector<Complex>* inValue   = NULL;

  for(size_t step = 0; step < maxIt; step++){
    // Formulations //
    wave = new FormulationSteadyWaveScalar<Complex>(*domain, k, order);

    // System //
    system = new System<Complex>(*wave);

    // Init DDM Border Dofs
    if(step == 0){
      solution = new map<Dof, Complex>;
      ddmG     = new map<Dof, Complex>;

      initMap(*system, *ddmBorder, *solution);
      initMap(*system, *ddmBorder, *ddmG);
    }

    // Init MPI Buffers
    if(step == 0){
      outEntity = new vector<int>(ddmG->size(), 0);
      outType   = new vector<int>(ddmG->size(), 0);
      outValue  = new vector<Complex>(ddmG->size(), 0);

      inEntity = new vector<int>(ddmG->size(), 0);
      inType   = new vector<int>(ddmG->size(), 0);
      inValue  = new vector<Complex>(ddmG->size(), 0);
    }

    // Constraint
    if(myId == 0)
      SystemHelper<Complex>::dirichlet(*system, *source, fSource);

    // Assemble //
    // Volume
    system->assemble();

    // Neumann terms
    neumann = new FormulationNeumann(*infinity, k, order);
    system->addBorderTerm(*neumann);

    // DDM terms
    if(ddmType == emdaType)
      ddm = new FormulationEMDA(*ddmBorder, k, chi, order, *ddmG);
    else if(ddmType == oo2Type)
      ddm = new FormulationOO2(*ddmBorder, ooA, ooB, order, *ddmG);
    else
      throw Exception("Unknown %s DDM border term", ddmType.c_str());

    system->addBorderTerm(*ddm);

    // Solve //
    system->solve();

    // Get DDM Border Solution //
    system->getSolution(*solution, 0);

    try{
      if(option.getValue("-disp").size() > 1)
        displaySolution(*solution, atoi(option.getValue("-disp")[1].c_str()),
                        step, "u");
    }
    catch(...){
    }

    // Update G //
    const FunctionSpaceScalar& ddmFSpace =
      static_cast<const FunctionSpaceScalar&>(ddm->fs());

    if(ddmType == emdaType)
      upDdm  = new FormulationUpdateEMDA(ddmFSpace, k, chi, *solution, *ddmG);
    else if(ddmType == oo2Type)
      upDdm  = new FormulationUpdateOO2(ddmFSpace, ooA, ooB, *solution, *ddmG);
    else
      throw Exception("Unknown %s DDM border term", ddmType.c_str());

    update = new System<Complex>(*upDdm);

    update->assemble();
    update->solve();
    update->getSolution(*ddmG, 0);

    // Serialize ddmG
    serialize(*ddmG, *outEntity, *outType, *outValue);

    // Exchange
    exchange(myId,
             *outEntity, *outType, *outValue, *inEntity, *inType, *inValue);

    // ddmG is new ddmG : unserialize
    unserialize(*ddmG, *inEntity, *inType, *inValue);

    // Write Solution //
    if(step == maxIt - 1){
      stringstream feSolName;
      FEMSolution<Complex> feSol;
      system->getSolution(feSol);

      feSolName << "ddm" << myId;
      feSol.write(feSolName.str());
    }

    // Clean //
    delete update;
    delete upDdm;

    delete ddm;
    if(neumann)
      delete neumann;

    delete wave;
    delete system;
  }

  // Finalize //
  delete domain;
  if(source)
    delete source;
  delete infinity;
  delete ddmBorder;

  delete solution;
  delete ddmG;

  delete outEntity;
  delete outType;
  delete outValue;

  delete inEntity;
  delete inType;
  delete inValue;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-max,-ddm,-chi,-lc,-disp");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
