#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "WriterMsh.h"

using namespace std;

void WriterMsh::writeNodalValuesHeader(const string name) const{
  *out << "$ElementNodeData"   << endl;

  if(isNodal)
    *out << "1"                        << endl  // 1 string tag
         << "\"" << name << "\""       << endl; // (name)

  else
    *out << "2"                        << endl  // 2 string tag
         << "\"" << name << "\""       << endl  // (name)
         << "\"interpolation scheme\"" << endl; // (interpolation scheme)

  *out << "1"                          << endl  // 1 real tag
       << "0"                          << endl  // (time value)
       << "3"                          << endl  // 3 integer tag
       << "0"                          << endl; // (time step index)

  if(isScalar)
    *out << "1" << endl;                // (number of field -- scalar)
  else
    *out << "3" << endl;                // (number of field -- vector)

  *out << E << endl;                    // (number of element)
}

void WriterMsh::writeNodalValuesFromNode(void) const{
  for(size_t i = 0; i < E; i++){
    *out << (*element)[i]->getNum()         << " "
         << (*element)[i]->getNumVertices() << " ";

    const size_t M      = (*element)[i]->getNumVertices();
    MElement* myElement = const_cast<MElement*>((*element)[i]);

    for(size_t j = 0; j < M; j++){
      const size_t id = myElement->getVertex(j)->getNum() - 1;
      // Note: getNum() ranges from *1* to MAX
      //   --> we need to substract 1 !!

      if(isScalar)
        *out << (*scalarValue)[id] << " ";
      else
        *out << (*vectorValue)[id](0) << " "
             << (*vectorValue)[id](1) << " "
             << (*vectorValue)[id](2) << " ";
    }

    *out << endl;
  }
}

void WriterMsh::writeNodalValuesFromSol(void) const{
  // Lagrange Basis Size //
  const size_t nCoef = lBasis->getNFunction();

  // Scalar FS ? //
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;

  // Temporary //
  size_t globalId;

  if(isScalar)
    fsScalar = static_cast<const FunctionSpaceScalar*>(fs);

  else
    fsVector = static_cast<const FunctionSpaceVector*>(fs);

  // Iterate on Element //
  for(size_t i = 0; i < E; i++){
    *out << (*element)[i]->getNum() << " "
         << nCoef                   << " ";

    // Get Element GoD
    const GroupOfDof& god = fs->getGoDFromElement(*(*element)[i]);

    // Get Dof
    const vector<const Dof*>& dof  = god.getDof();
    const size_t              size = dof.size();

    // Get Coef In FS Basis
    vector<double> fsCoef(size);
    for(size_t j = 0; j < size; j++){
      // Dof Global ID
      globalId = dofM->getGlobalId(*dof[j]);

      // If non fixed Dof: look in Solution
      if(globalId != DofManager::isFixedId())
        fsCoef[j] = (*sol)(globalId);

      // If Dof is fixed: get fixed value
      else
        fsCoef[j] = dofM->getValue(*dof[j]);
    }

    // Get Coef In Lagrange Basis
    if(isScalar){
      vector<double> lCoef =
        lBasis->project(*(*element)[i], fsCoef, *fsScalar);

      for(size_t j = 0; j < nCoef; j++)
        *out << lCoef[j] << " ";
    }

    else{
      vector<fullVector<double> > lCoef =
        lBasis->project(*(*element)[i], fsCoef, *fsVector);

      for(size_t j = 0; j < nCoef; j++)
        *out << lCoef[j](0) << " "
             << lCoef[j](1) << " "
             << lCoef[j](2) << " ";
    }

    *out << endl;
  }
}

void WriterMsh::writeNodalValuesFooter(void) const{
  *out << "$EndElementNodeData" << endl;
}
