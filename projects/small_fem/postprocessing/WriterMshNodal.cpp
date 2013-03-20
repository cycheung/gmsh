#include "DofFixedException.h"
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
  for(unsigned int i = 0; i < E; i++){
    *out << (*element)[i]->getNum()         << " "
	 << (*element)[i]->getNumVertices() << " ";

    const unsigned int M = (*element)[i]->getNumVertices();
    MElement* myElement = const_cast<MElement*>((*element)[i]);

    for(unsigned int j = 0; j < M; j++){
      const unsigned int id = myElement->getVertex(j)->getNum() - 1;
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
  const unsigned int nCoef = lBasis->getNFunction();

  // Scalar FS ? //
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;

  if(isScalar)
    fsScalar = static_cast<const FunctionSpaceScalar*>(fs);

  else
    fsVector = static_cast<const FunctionSpaceVector*>(fs);

  // Iterate on Element //
  for(unsigned int i = 0; i < E; i++){
    *out << (*element)[i]->getNum() << " "
	 << nCoef                   << " ";

    // Get Element GoD
    const GroupOfDof& god = fs->getGoDFromElement(*(*element)[i]);

    // Get Dof
    const vector<const Dof*>& dof  = god.getAll();
    const unsigned int        size = dof.size();

    // Get Coef In FS Basis
    vector<double> fsCoef(size);
    for(unsigned int j = 0; j < size; j++)
          try{
            // Look in Solution
            fsCoef[j] = (*sol)(dofM->getGlobalId(*dof[j]));
          }

          catch(DofFixedException& fixedDof){
            // If Dos is fixed, look in 'fixedDof'
            fsCoef[j] = fixedDof.getValue();
          }

    // Get Coef In Lagrange Basis
    if(isScalar){
      vector<double> lCoef =
	lBasis->project(*(*element)[i], fsCoef, *fsScalar);

      for(unsigned int j = 0; j < nCoef; j++)
	*out << lCoef[j] << " ";
    }

    else{
      vector<fullVector<double> > lCoef =
	lBasis->project(*(*element)[i], fsCoef, *fsVector);

      for(unsigned int j = 0; j < nCoef; j++)
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
