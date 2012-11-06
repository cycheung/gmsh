#include <sstream>

#include "Exception.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "LagrangeGenerator.h"
#include "WriterMsh.h"

using namespace std;

WriterMsh::WriterMsh(void){
}

WriterMsh::~WriterMsh(void){
}

void WriterMsh::write(const std::string name) const{
  stringstream fileName; 
  fileName << name << ".msh";
  
  out = new ofstream;
  out->open(fileName.str().c_str());
  
  if(!hasDomain){
    *out << "No Domain has been given !" << endl;
  }
  
  else{
    writeHeader();
    writeNodes();
    writeElements();
    
    if(!isNodal){
      lBasis = LagrangeGenerator::generate((*element)[0]->getType(),
					   fs->getOrder());
      writeInterpolationScheme();
    }

    if(hasValue){
      writeNodalValuesHeader(name);
      
      if(isNodal)
	writeNodalValuesFromNode();

      else
	writeNodalValuesFromSol();
  
      writeNodalValuesFooter();  
    }
  }

  out->close();
  delete out;

  if(!isNodal)
    delete lBasis;
}

void WriterMsh::writeHeader(void) const{
  *out << "$MeshFormat" << endl
       << "2.2 0 8" << endl
       << "$EndMeshFormat" << endl; 
}

void WriterMsh::writeNodes(void) const{
  *out << "$Nodes" << endl
       << N << endl;

  for(int i = 0; i < N; i++){
    *out << (*node)[i]->getNum() << " "
	 << (*node)[i]->x()      << " "
	 << (*node)[i]->y()      << " "
	 << (*node)[i]->z()      << endl;
  }

  *out << "$EndNodes" << endl;
}

void WriterMsh::writeElements(void) const{
  *out << "$Elements" << endl
       << E << endl;
  
  for(int i = 0; i < E; i++){
    *out << (*element)[i]->getNum()        << " " 
	 << (*element)[i]->getTypeForMSH() 
	 << " 2 1 1" << " "; 
           // 2 Tags --> (1 physical entity, 1 elementary geometry) 

    const int M = (*element)[i]->getNumVertices();
    MElement* myElement = 
      const_cast<MElement*>((*element)[i]);

    for(int j = 0; j < M; j++)
      *out << myElement->getVertex(j)->getNum() << " ";
    
    *out << endl;
  }
  
  *out << "$EndElements" << endl;
}

void WriterMsh::writeInterpolationScheme(void) const{ 
  // Up to now, adaptive view with only scalar values
  if(!isScalar)
    throw Exception("Adaptive View with Scalar Only");

  // Some Temp Value
  const fullMatrix<double>& coef = lBasis->getCoefficient();
  const fullMatrix<double>& mono = lBasis->getMonomial();

  const unsigned int nRowCoef = coef.size1();
  const unsigned int nColCoef = coef.size2();
  
  const unsigned int nRowMono = mono.size1();
  const unsigned int nColMono = mono.size2();

  // Up to now, we suppose *ONE* topology
  *out << "$InterpolationScheme"     << endl
       << "\"interpolation scheme\"" << endl
       << "1"                        << endl 
       << (*element)[0]->getType()   << endl;
  
  if(isScalar)
    *out << "2" << endl; // 2 Matrices for Scalar
  
  else
    *out << "6" << endl; // 6 ?? for Vector
  
  if(isScalar){
    // Scalar Coef Matrix
    *out << nRowCoef << " "
	 << nColCoef << endl;

    for(unsigned int i = 0; i < nRowCoef; i++){
      for(unsigned int j = 0; j < nColCoef; j++){
	*out << coef(i, j);
	
	if(j < nColCoef - 1)
	  *out << " ";
      
	else
	  *out << endl;
      }
    }

    // Scalar Mono Matrix
    *out << nRowMono << " "
	 << nColMono << endl;

    for(unsigned int i = 0; i < nRowMono; i++){
      for(unsigned int j = 0; j < nColMono; j++){
	*out << mono(i, j);
      
	if(j < nColMono - 1)
	  *out << " ";
	
	else
	  *out << endl;
      }
    }
  }

  else{
  }

  *out << "$EndInterpolationScheme" << endl;
}

void WriterMsh::writeNodalValuesHeader(const string name) const{
  *out << "$ElementNodeData"   << endl;

  if(isNodal)
    *out << "1"                  << endl  // 1 string tag
	 << "\"" << name << "\"" << endl; // (name)

  else
    *out << "2"                        << endl  // 2 string tag
	 << "\"" << name << "\""       << endl  // (name)
	 << "\"interpolation scheme\"" << endl; // (interpolation scheme)
    
  *out << "1" << endl  // 1 real tag 
       << "0" << endl  // (time value)
       << "3" << endl  // 3 integer tag
       << "0" << endl; // (time step index)
    
  if(isScalar)
    *out << "1" << endl;                // (number of field -- scalar)
  else
    *out << "3" << endl;                // (number of field -- vector)
  
  *out << E << endl;                    // (number of element)
}

void WriterMsh::writeNodalValuesFromNode(void) const{
  for(int i = 0; i < E; i++){
    *out << (*element)[i]->getNum()         << " " 
	 << (*element)[i]->getNumVertices() << " ";
    
    const int M = (*element)[i]->getNumVertices();
    MElement* myElement = 
      const_cast<MElement*>((*element)[i]);

    for(int j = 0; j < M; j++){
      const int id = myElement->getVertex(j)->getNum() - 1;
      // Note: getNum() ranges from *1* to MAX
      //   --> we need to substract 1 !!

      if(isScalar)
	*out << (*nodalScalarValue)[id] << " ";
      else
	*out << (*nodalVectorValue)[id](0) << " "
	     << (*nodalVectorValue)[id](1) << " "
	     << (*nodalVectorValue)[id](2) << " ";
    }
    
    *out << endl;
  }
}

void WriterMsh::writeNodalValuesFromSol(void) const{
  // Lagrange Basis Size //
  const unsigned int nCoef = lBasis->getSize();

  // Scalar FS ? //
  const FunctionSpaceScalar* fsScalar = NULL;
  const FunctionSpaceVector* fsVector = NULL;
 
  if(isScalar)
    fsScalar = static_cast<const FunctionSpaceScalar*>(fs);

  else
    fsVector = static_cast<const FunctionSpaceVector*>(fs);

  // Iterate on Element //
  for(int i = 0; i < E; i++){
    *out << (*element)[i]->getNum() << " " 
	 << nCoef                   << " ";
    
    // Get Element GoD
    const GroupOfDof& god = fs->getGoDFromElement(*(*element)[i]);

    // Get Dof
    const vector<const Dof*>& dof  = god.getAll();
    const unsigned int        size = dof.size();
    
    // Get Coef In FS Basis
    fullVector<double> fsCoef(size);
    for(unsigned int j = 0; j < size; j++)
      // Look in Solution
      fsCoef(j) = 
	(*sol)(dofM->getGlobalId(*dof[j])); 

    // Get Coef In Lagrange Basis
    if(isScalar){
      fullVector<double> lCoef = 
	lBasis->project(fsCoef, fsScalar->getLocalFunctions(*(*element)[i]));
      
      for(unsigned int j = 0; j < nCoef; j++)
	*out << lCoef(j) << " ";
    
    }

    *out << endl;
  }
}

void WriterMsh::writeNodalValuesFooter(void) const{
  *out << "$EndElementNodeData" << endl;
}
