#include "FEMSolution.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "BasisLagrange.h"
#include "BasisGenerator.h"

#include "Exception.h"

using namespace std;

FEMSolution::FEMSolution(void){
  pView = new PViewDataGModel(PViewDataGModel::ElementNodeData);
}

FEMSolution::~FEMSolution(void){
  pView->destroyData();
  delete pView;
}

void FEMSolution::addCoefficients(size_t step,
                                  double time,
                                  const FunctionSpace& fs,
                                  const DofManager& dofM,
                                  const fullVector<double>& coef){
  // Get Support and GModel //
  const vector<const MElement*>& element = fs.getSupport().getAll();
  const size_t                  nElement = element.size();
  GModel&                          model = fs.getSupport().getMesh().getModel();


  // Lagrange Basis & Interpolation matrices //
  BasisLagrange* lagrange = static_cast<BasisLagrange*>
    (BasisGenerator::generate(element[0]->getType(),
                              0,
                              fs.getBasis(0).getOrder(),
                              "lagrange"));

  pView->setInterpolationMatrices(element[0]->getType(),
                                  lagrange->getCoefficient(),
                                  lagrange->getMonomial());

  // Map with (Element Id, Lagrange coefficients) //
  map<int, vector<double> > data;

  // Scalar of Vectorial Field ? //
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;
  size_t nComp;

  if(fs.isScalar()){
    fsScalar = static_cast<const FunctionSpaceScalar*>(&fs);
    nComp = 1;
  }

  else{
    fsVector = static_cast<const FunctionSpaceVector*>(&fs);
    nComp = 3;
  }

  // Iterate on Element //
  for(size_t i = 0; i < nElement; i++){
    // Get Element GoD
    const GroupOfDof& god = fs.getGoDFromElement(*element[i]);

    // Get Dof
    const vector<Dof>& dof  = god.getDof();
    const size_t       size = dof.size();

    // Get Coef In FS Basis
    vector<double> fsCoef(size);
    for(size_t j = 0; j < size; j++){
      // Dof Global ID
      size_t globalId = dofM.getGlobalId(dof[j]);

      // If non fixed Dof: look in Solution
      if(globalId != DofManager::isFixedId())
        fsCoef[j] = coef(globalId);

      // If Dof is fixed: get fixed value
      else
        fsCoef[j] = dofM.getValue(dof[j]);
    }

    // Get Coef In Lagrange Basis
    vector<double> lCoef;
    if(fsScalar)
      lCoef = lagrange->project(*element[i], fsCoef, *fsScalar);

    else
      lCoef = lagrange->project(*element[i], fsCoef, *fsVector);

    // Add in map
    data.insert(pair<int, vector<double> >(element[i]->getNum(), lCoef));
  }

  // Add map to PView //
  pView->addData(&model, data, step, time, 0, nComp);

  // Clean //
  delete lagrange;
}

void FEMSolution::addCoefficients(size_t step,
                                  double time,
                                  const FunctionSpace& fs,
                                  const DofManager& dofM,
                                  const fullVector<std::complex<double> >& coef)
{
  // Complex coefficients lead to two map (one real, one imaginary) //
  // So steps goes in pair (even: real part; odd: imaginary part)   //

  // Check if step is even //
  if(step % 2 != 0)
    throw Exception("%s: %s -- %s, %s", "FEMSolution",
                    "cannot add coefficients",
                    "with complex values two 'steps' are inserted",
                    "so 'step' must be an even number");

  // Split real and imaginary parts //
  const size_t size = coef.size();
  fullVector<double> real(size);
  fullVector<double> imag(size);

  for(size_t i = 0; i < size; i++)
    real(i) = coef(i).real();

  for(size_t i = 0; i < size; i++)
    imag(i) = coef(i).imag();

  // Add real & imaginary parts //
  addCoefficients(step + 0, time, fs, dofM, real);
  addCoefficients(step + 1, time, fs, dofM, imag);
}

void FEMSolution::clear(void){
  pView->destroyData();
}

void FEMSolution::write(std::string fileName) const{
  pView->setName(fileName);
  pView->writeMSH(fileName + ".msh");
}
