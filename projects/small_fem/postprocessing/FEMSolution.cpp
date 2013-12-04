#include "FEMSolution.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "BasisLagrange.h"
#include "BasisGenerator.h"

#include "Exception.h"

using namespace std;

template<>
void FEMSolution<double>::
addCoefficients(size_t step,
                double time,
                const FunctionSpace& fs,
                const DofManager<double>& dofM,
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
      if(globalId != DofManager<double>::isFixedId())
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

template<>
void FEMSolution<complex<double> >::
addCoefficients(size_t step,
                double time,
                const FunctionSpace& fs,
                const DofManager<complex<double> >& dofM,
                const fullVector<complex<double> >& coef){

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
  // Real And Imaginary parts
  map<int, vector<double> > real;
  map<int, vector<double> > imag;

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
    vector<complex<double> > fsCoef(size);
    for(size_t j = 0; j < size; j++){
      // Dof Global ID
      size_t globalId = dofM.getGlobalId(dof[j]);

      // If non fixed Dof: look in Solution
      if(globalId != DofManager<complex<double> >::isFixedId())
        fsCoef[j] = coef(globalId);

      // If Dof is fixed: get fixed value
      else
        fsCoef[j] = dofM.getValue(dof[j]);
    }

    // Split fsCoef
    vector<double> fsCoefReal(size);
    vector<double> fsCoefImag(size);

    for(size_t j = 0; j < size; j++)
      fsCoefReal[j] = fsCoef[j].real();

    for(size_t j = 0; j < size; j++)
      fsCoefImag[j] = fsCoef[j].imag();

    // Get Coef In Lagrange Basis
    vector<double> lCoefReal;
    vector<double> lCoefImag;

    if(fsScalar){
      lCoefReal = lagrange->project(*element[i], fsCoefReal, *fsScalar);
      lCoefImag = lagrange->project(*element[i], fsCoefImag, *fsScalar);
    }

    else{
      lCoefReal = lagrange->project(*element[i], fsCoefReal, *fsVector);
      lCoefImag = lagrange->project(*element[i], fsCoefImag, *fsVector);
    }

    // Add in map
    real.insert(pair<int, vector<double> >(element[i]->getNum(), lCoefReal));
    imag.insert(pair<int, vector<double> >(element[i]->getNum(), lCoefImag));
  }

  // Add map to PView //
  pView->addData(&model, real, 2 * step + 0, time, 0, nComp);
  pView->addData(&model, imag, 2 * step + 1, time, 0, nComp);

  // Clean //
  delete lagrange;
}
