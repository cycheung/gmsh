////////////////////////////////////////////////
// Templates Implementations for FEMSolution: //
// Inclusion compilation model                //
//                                            //
// Damn you gcc: we want 'export' !           //
////////////////////////////////////////////////

template<typename scalar>
FEMSolution<scalar>::FEMSolution(void){
  pView = new PViewDataGModel(PViewDataGModel::ElementNodeData);
}

template<typename scalar>
FEMSolution<scalar>::~FEMSolution(void){
  pView->destroyData();
  delete pView;
}

template<typename scalar>
void FEMSolution<scalar>::clear(void){
  pView->destroyData();
}

template<typename scalar>
void FEMSolution<scalar>::write(std::string fileName) const{
  pView->setName(fileName);
  pView->writeMSH(fileName + ".msh");
}
