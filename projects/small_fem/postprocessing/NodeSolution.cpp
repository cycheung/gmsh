#include "NodeSolution.h"

#include "Exception.h"

using namespace std;

NodeSolution::NodeSolution(void){
  pView = new PViewDataGModel(PViewDataGModel::NodeData);
}

NodeSolution::~NodeSolution(void){
  pView->destroyData();
  delete pView;
}

void NodeSolution::write(std::string fileName) const{
  pView->setName(fileName);
  pView->writeMSH(fileName + ".msh");
}


void NodeSolution::
addNodeValue(size_t step,
             double time,
             const Mesh& mesh,
             std::map<const MVertex*, std::complex<double> >& data){

  // GModel //
  GModel& model = mesh.getModel();

  // Map with (Vertex Id, Node Value) //
  map<int, vector<double> > gmshDataReal;
  map<int, vector<double> > gmshDataImag;

  // Scalar of Vectorial Field ? //
  size_t nComp;
  nComp = 1; // Up to now: scalar field

  // Populate gmshData //
  map<const MVertex*, complex<double> >::iterator it  = data.begin();
  map<const MVertex*, complex<double> >::iterator end = data.end();

  vector<double> tmpReal(1);
  vector<double> tmpImag(1);

  for(; it != end; it++){
    tmpReal[0] = it->second.real();
    tmpImag[0] = it->second.imag();

    gmshDataReal.insert
      (pair<int, vector<double> >(it->first->getNum(), tmpReal));
    gmshDataImag.insert
      (pair<int, vector<double> >(it->first->getNum(), tmpImag));
  }

  // Add map to PView //
  pView->addData(&model, gmshDataReal, 2 * step + 0, time, 0, nComp);
  pView->addData(&model, gmshDataImag, 2 * step + 1, time, 0, nComp);
}
