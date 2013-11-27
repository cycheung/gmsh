#include "ElementSolution.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "BasisLagrange.h"
#include "BasisGenerator.h"

#include "Exception.h"

using namespace std;

ElementSolution::ElementSolution(void){
  pView = new PViewDataGModel(PViewDataGModel::ElementData);
}

ElementSolution::~ElementSolution(void){
  pView->destroyData();
  delete pView;
}

void ElementSolution::clear(void){
  pView->destroyData();
}

void ElementSolution::addValues(size_t step,
                                double time,
                                const GroupOfElement& goe,
                                const std::vector<double>& value){
  // Get Elements and GModel //
  const vector<const MElement*>& element = goe.getAll();
  GModel&                          model = goe.getMesh().getModel();

  // Map all Element Ids with its value //
  const size_t nElement = element.size();
  vector<double>            tmp(1);
  map<int, vector<double> > data;

  for(size_t i = 0; i < nElement; i++){
    tmp[0] = value[i];
    data.insert(pair<int, vector<double> >(element[i]->getNum(), tmp));
  }

  // Add map to PView //
  pView->addData(&model, data, step, time, 0, 1);
}

void ElementSolution::write(std::string fileName) const{
  pView->setName(fileName);
  pView->writeMSH(fileName + ".msh");
}
