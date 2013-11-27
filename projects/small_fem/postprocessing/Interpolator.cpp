#include "Interpolator.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

using namespace std;

Interpolator::Interpolator(void){
}

Interpolator::~Interpolator(void){
}

void Interpolator::interpolate(const FunctionSpace& fs,
                               const DofManager& dofM,
                               const fullVector<double>& coef,
                               const fullMatrix<double>& point,
                               fullMatrix<double>& values){
  // Get GModel //
  GModel&    model = fs.getSupport().getMesh().getModel();
  const size_t dim = model.getDim();

  // Scalar or Vector ?
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;
  const bool                 isScalar = fs.isScalar();

  if(isScalar)
    fsScalar = static_cast<const FunctionSpaceScalar*>(&fs);
  else
    fsVector = static_cast<const FunctionSpaceVector*>(&fs);

  // Alloc values //
  const size_t nPoint = point.size1();

  if(isScalar)
    values.resize(nPoint, 1);
  else
    values.resize(nPoint, 3);

  // Iterate on 'point'
  for(size_t i = 0; i < nPoint; i++){
    // Search element containg this point
    SPoint3   thisPoint(point(i, 0), point(i, 1), point(i, 2));
    MElement* element = model.getMeshElementByCoord(thisPoint, dim, true);

    // WARNING: if no element found, set 'values' to zero
    if(!element){
      values(i, 0) = 0;

      if(!isScalar){
        values(i, 1) = 0;
        values(i, 2) = 0;
      }
    }

    else{
      // Get GroupOfDof related to this Element
      const GroupOfDof& god = fs.getGoDFromElement(*element);

      // Get Dof
      const vector<Dof>& dof  = god.getDof();
      const size_t       size = dof.size();

      // Get Coef
      vector<double> thisCoef(size);
      for(size_t k = 0; k < size; k++){
        // Dof Global ID
        const size_t globalId = dofM.getGlobalId(dof[k]);

        // If non fixed Dof: look in Solution
        if(globalId != DofManager::isFixedId())
          thisCoef[k] = coef(globalId);

        // If Dof is fixed: get fixed value
        else
          thisCoef[k] = dofM.getValue(dof[k]);
      }

      // Get Node coordinate
      fullVector<double> xyz(3);
      xyz(0) = point(i, 0);
      xyz(1) = point(i, 1);
      xyz(2) = point(i, 2);

      // Interpolate (AT LAST !!)
      if(isScalar){
        values(i, 0) = fsScalar->interpolate(*element, thisCoef, xyz);
      }

      else{
        fullVector<double> tmp = fsVector->interpolate(*element, thisCoef, xyz);

        values(i, 0) = tmp(0);
        values(i, 1) = tmp(1);
        values(i, 2) = tmp(2);
      }
    }
  }
}
