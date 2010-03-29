
// Gmsh - Copyright (C) 1997-2010 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.
//
// Contributed by Matti Pellikka <matti.pellikka@tut.fi>.

#include "ChainComplex.h"

#if defined(HAVE_KBIPACK)

ChainComplex::ChainComplex(CellComplex* cellComplex, int domain)
{ 
  _dim = cellComplex->getDim();
  _cellComplex = cellComplex;
  
  for(int i = 0; i < 5; i++){
    _HMatrix[i] = NULL;
    _kerH[i] = NULL;
    _codH[i] = NULL;
    _JMatrix[i] = NULL;
    _QMatrix[i] = NULL;
    _Hbasis[i] = NULL;
  }

  int lastCols = 0;
  for(int dim = 0; dim < 4; dim++){
    unsigned int cols = cellComplex->getSize(dim);
    unsigned int rows = 0;
    
    int index = 1;
    // ignore cells depending on domain
    for(CellComplex::citer cit = cellComplex->firstCell(dim); 
	cit != cellComplex->lastCell(dim); cit++){
      Cell* cell = *cit;
      cell->setIndex(0);
      cols--;
      if((domain == 0 && !cell->inSubdomain()) || domain == 1 
	 || (domain == 2 && cell->inSubdomain()) ){
        cols++;
	cell->setIndex(index);
	index++;
	_cellIndices[dim][cell] = cell->getIndex();
      }
    }
    if(dim > 0) rows = lastCols;
    lastCols = cols;
    
    if(cols == 0){ // no dim-cells, no map
      //_HMatrix[dim] = create_gmp_matrix_zero(rows, 1);
      _HMatrix[dim] = NULL;
    }
    else if(rows == 0){ // no dim-1-cells, maps everything to zero
      _HMatrix[dim] = create_gmp_matrix_zero(1, cols);
      //_HMatrix[dim] = NULL;
    }
    
    else{
      mpz_t elem;
      mpz_init(elem);
      _HMatrix[dim] = create_gmp_matrix_zero(rows, cols);
      for( std::set<Cell*, Less_Cell>::iterator cit = 
	     cellComplex->firstCell(dim);
	   cit != cellComplex->lastCell(dim); cit++){
        Cell* cell = *cit;
        if( (domain == 0 && !cell->inSubdomain()) || domain == 1 
	    || (domain == 2 && cell->inSubdomain()) ){
          for(Cell::biter it = cell->firstBoundary();
	      it != cell->lastBoundary(); it++){
            Cell* bdCell = (*it).first;
            if((domain == 0 && !bdCell->inSubdomain()) || domain == 1 
	       || (domain == 2 && cell->inSubdomain()) ){
              int old_elem = 0;

              if(bdCell->getIndex() > (int)gmp_matrix_rows( _HMatrix[dim]) 
		 || bdCell->getIndex() < 1 
                 || cell->getIndex() > (int)gmp_matrix_cols( _HMatrix[dim]) 
		 || cell->getIndex() < 1){
                printf("Warning: Index out of bound! HMatrix: %d. \n", dim);
              }
              else{
                gmp_matrix_get_elem(elem, bdCell->getIndex(), 
				    cell->getIndex(), _HMatrix[dim]);
                old_elem = mpz_get_si(elem);
                mpz_set_si(elem, old_elem + (*it).second);
                if( abs((old_elem + (*it).second)) > 1){
		  //printf("Incidence index: %d, in HMatrix: %d. \n", (old_elem + (*it).second), dim);
                }
                gmp_matrix_set_elem(elem, bdCell->getIndex(), 
				    cell->getIndex(), _HMatrix[dim]);
              }
            }
          }
        }
      } 
      mpz_clear(elem); 
    }

    _kerH[dim] = NULL;
    _codH[dim] = NULL;
    _JMatrix[dim] = NULL;
    _QMatrix[dim] = NULL;
    _Hbasis[dim] = NULL;     
  }
}

ChainComplex::~ChainComplex()
{
  for(int i = 0; i < 5; i++){
    destroy_gmp_matrix(_HMatrix[i]);
    destroy_gmp_matrix(_kerH[i]);
    destroy_gmp_matrix(_codH[i]);
    destroy_gmp_matrix(_JMatrix[i]);
    destroy_gmp_matrix(_QMatrix[i]);
    destroy_gmp_matrix(_Hbasis[i]);
  }
}

void ChainComplex::KerCod(int dim)
{ 
  if(dim < 0 || dim > 3 || _HMatrix[dim] == NULL) return;
  
  gmp_matrix* HMatrix 
    = copy_gmp_matrix(_HMatrix[dim], 1, 1, 
		      gmp_matrix_rows(_HMatrix[dim]),
		      gmp_matrix_cols(_HMatrix[dim]));
  
  gmp_normal_form* normalForm 
    = create_gmp_Hermite_normal_form(HMatrix, NOT_INVERTED, INVERTED);
  //printMatrix(normalForm->left);
  //printMatrix(normalForm->canonical);
  //printMatrix(normalForm->right);
  
  int minRowCol = std::min(gmp_matrix_rows(normalForm->canonical), 
			   gmp_matrix_cols(normalForm->canonical));
  int rank = 0;
  mpz_t elem;
  mpz_init(elem);
  
  // find the rank
  while(rank < minRowCol){
    gmp_matrix_get_elem(elem, rank+1, rank+1, normalForm->canonical);
    if(mpz_cmp_si(elem,0) == 0) break;
    rank++;
  }
  
  if(rank != (int)gmp_matrix_cols(normalForm->canonical)){
    _kerH[dim] 
      = copy_gmp_matrix(normalForm->right, 1, rank+1, 
			gmp_matrix_rows(normalForm->right),
			gmp_matrix_cols(normalForm->right));
  }
  
  if(rank > 0){
     _codH[dim] = 
       copy_gmp_matrix(normalForm->canonical, 1, 1,
		       gmp_matrix_rows(normalForm->canonical), rank);
     gmp_matrix_left_mult(normalForm->left, _codH[dim]);
  }
  
  mpz_clear(elem);
  destroy_gmp_normal_form(normalForm);
  
  return;
}

//j:B_k->Z_k
void ChainComplex::Inclusion(int lowDim, int highDim)
{
  if(getKerHMatrix(lowDim) == NULL 
     || getCodHMatrix(highDim) == NULL 
     || abs(lowDim-highDim) != 1) return;
  
  gmp_matrix* Zbasis = 
    copy_gmp_matrix(_kerH[lowDim], 1, 1,
		    gmp_matrix_rows(_kerH[lowDim]), 
		    gmp_matrix_cols(_kerH[lowDim]));
  gmp_matrix* Bbasis 
    = copy_gmp_matrix(_codH[highDim], 1, 1,
		      gmp_matrix_rows(_codH[highDim]), 
		      gmp_matrix_cols(_codH[highDim]));
  
  
  int rows = gmp_matrix_rows(Bbasis);
  int cols = gmp_matrix_cols(Bbasis);
  if(rows < cols) return;
  
  rows = gmp_matrix_rows(Zbasis);
  cols = gmp_matrix_cols(Zbasis);
  if(rows < cols) return;
  
  // inv(U)*A*inv(V) = S
  gmp_normal_form* normalForm 
    = create_gmp_Smith_normal_form(Zbasis, INVERTED, INVERTED);
  
  mpz_t elem;
  mpz_init(elem);
  
  for(int i = 1; i <= cols; i++){
  
    gmp_matrix_get_elem(elem, i, i, normalForm->canonical);
    if(mpz_cmp_si(elem,0) == 0){
      destroy_gmp_normal_form(normalForm);
      return;
    }
  }
  
  gmp_matrix_left_mult(normalForm->left, Bbasis); 
  
  gmp_matrix* LB = copy_gmp_matrix(Bbasis, 1, 1, 
				   gmp_matrix_cols(Zbasis), 
				   gmp_matrix_cols(Bbasis));
  destroy_gmp_matrix(Bbasis);
  
  rows = gmp_matrix_rows(LB);
  cols = gmp_matrix_cols(LB);
  
  mpz_t divisor;
  mpz_init(divisor);
  mpz_t remainder;
  mpz_init(remainder);
  mpz_t result;
  mpz_init(result);
  
  for(int i = 1; i <= rows; i++){
    gmp_matrix_get_elem(divisor, i, i, normalForm->canonical);
    for(int j = 1; j <= cols; j++){
      gmp_matrix_get_elem(elem, i, j, LB);
      mpz_cdiv_qr(result, remainder, elem, divisor);
      if(mpz_cmp_si(remainder, 0) == 0){
        gmp_matrix_set_elem(result, i, j, LB);
      }
      else return;
    }
  }
  
  gmp_matrix_left_mult(normalForm->right, LB);
  
  setJMatrix(lowDim, LB);
  
  mpz_clear(elem);
  mpz_clear(divisor);
  mpz_clear(result);
  destroy_gmp_normal_form(normalForm);
}

void ChainComplex::Quotient(int dim)
{
  if(dim < 0 || dim > 4 || _JMatrix[dim] == NULL) return;
  
  gmp_matrix* JMatrix = 
    copy_gmp_matrix(_JMatrix[dim], 1, 1,
		    gmp_matrix_rows(_JMatrix[dim]), 
		    gmp_matrix_cols(_JMatrix[dim]));
  int rows = gmp_matrix_rows(JMatrix);
  int cols = gmp_matrix_cols(JMatrix);
  
  gmp_normal_form* normalForm = 
    create_gmp_Smith_normal_form(JMatrix, NOT_INVERTED, NOT_INVERTED);

  //printMatrix(normalForm->left);
  //printMatrix(normalForm->canonical);
  //printMatrix(normalForm->right);  
  
  mpz_t elem;
  mpz_init(elem);
    
  for(int i = 1; i <= cols; i++){
    gmp_matrix_get_elem(elem, i, i, normalForm->canonical);
    if(mpz_cmp_si(elem,0) == 0){
      destroy_gmp_normal_form(normalForm);
      return;
    }
    if(mpz_cmp_si(elem,1) > 0) _torsion[dim].push_back(mpz_get_si(elem));
  }
  
  int rank = cols - _torsion[dim].size();
  if(rows - rank > 0){
    gmp_matrix* Hbasis = 
      copy_gmp_matrix(normalForm->left, 1, rank+1, rows, rows);
    _QMatrix[dim] = Hbasis;
  }
  
  mpz_clear(elem);
  destroy_gmp_normal_form(normalForm);
  return; 
}

void ChainComplex::computeHomology(bool dual)
{  
  int lowDim = 0;
  int highDim = 0;
  int setDim = 0; 
  
  for(int i=-1; i < 4; i++){
    
    if(dual){
      lowDim = getDim()+1-i; 
      highDim = getDim()+1-(i+1);
      setDim = highDim;
      //KerCod(lowDim);
    }
    else{
      lowDim = i;
      highDim = i+1;
      setDim = lowDim;
      //KerCod(highDim);
    }
    
    //printf("Homology computation process: step %d of 4 \n", i+1);
    
    KerCod(highDim);
    
    // 1) no edges, but zero cells
    if(lowDim == 0 && !dual 
       &&  gmp_matrix_cols(getHMatrix(lowDim)) > 0 
       && getHMatrix(highDim) == NULL) {
      setHbasis( setDim, 
		 create_gmp_matrix_identity(gmp_matrix_cols(getHMatrix(lowDim))) );
    }
    else if(highDim == 0 && dual 
	    &&  gmp_matrix_rows(getHMatrix(highDim)) > 0 
	    && getHMatrix(lowDim) == NULL) {
      setHbasis( setDim, 
		 create_gmp_matrix_identity(gmp_matrix_rows(getHMatrix(highDim))) );
    }
    
    // 2) this dimension is empty
    else if(getHMatrix(lowDim) == NULL && getHMatrix(highDim) == NULL){
      setHbasis(setDim, NULL);
    }
    // 3) No higher dimension cells -> none of the cycles are boundaries
    else if(getHMatrix(highDim) == NULL){
      setHbasis( setDim, 
		 copy_gmp_matrix(getKerHMatrix(lowDim), 1, 1,
				 gmp_matrix_rows(getKerHMatrix(lowDim)), 
				 gmp_matrix_cols(getKerHMatrix(lowDim))) );
    }
    
   
    // 5) General case:
    //   1) Find the bases of boundaries B and cycles Z 
    //   2) find j: B -> Z and
    //   3) find quotient Z/j(B) 
    else {
      
      // 4) No lower dimension cells -> all chains are cycles
      if(getHMatrix(lowDim) == NULL){
        setKerHMatrix(lowDim, 
		      create_gmp_matrix_identity(gmp_matrix_rows(getHMatrix(highDim))) );
      }
      Inclusion(lowDim, highDim);
      Quotient(lowDim);
      
      if(getCodHMatrix(highDim) == NULL){
        setHbasis(setDim, 
		  copy_gmp_matrix(getKerHMatrix(lowDim), 1, 1,
				  gmp_matrix_rows(getKerHMatrix(lowDim)), 
				  gmp_matrix_cols(getKerHMatrix(lowDim))) );
      }  
      else if(getJMatrix(lowDim) == NULL || getQMatrix(lowDim) == NULL){
        setHbasis(setDim, NULL);
      } 
      else{
        setHbasis(setDim, 
		  copy_gmp_matrix(getKerHMatrix(lowDim), 1, 1, 
				  gmp_matrix_rows(getKerHMatrix(lowDim)), 
				  gmp_matrix_cols(getKerHMatrix(lowDim))) );
        
        gmp_matrix_right_mult(getHbasis(setDim), getQMatrix(lowDim));
      } 
    } 

    //destroy_gmp_matrix(getKerHMatrix(lowDim));
    //destroy_gmp_matrix(getCodHMatrix(lowDim));
    destroy_gmp_matrix(getJMatrix(lowDim));
    destroy_gmp_matrix(getQMatrix(lowDim));
    
    //setKerHMatrix(lowDim, NULL);
    //setCodHMatrix(lowDim, NULL);
    setJMatrix(lowDim, NULL);
    setQMatrix(lowDim, NULL);  
  } 
  return;
}


void ChainComplex::matrixTest()
{  
  const int rows = 3;
  const int cols = 6;
  
  long int elems[rows*cols];
  for(int i = 1; i<=rows*cols; i++) elems[i-1] = i;
  
  gmp_matrix* matrix = create_gmp_matrix_int(rows, cols, elems);
  
  gmp_matrix* copymatrix = copy_gmp_matrix(matrix, 3, 2, 3, 5);
  
  printMatrix(matrix);
  printMatrix(copymatrix);
}

std::vector<int> ChainComplex::getCoeffVector(int dim, int chainNumber)
{  
  std::vector<int> coeffVector;
  
  if(dim < 0 || dim > 4) return coeffVector;
  if(_Hbasis[dim] == NULL 
     || (int)gmp_matrix_cols(_Hbasis[dim]) < chainNumber) return coeffVector;
  
  int rows = gmp_matrix_rows(_Hbasis[dim]);
  
  int elemi;
  long int elemli;
  mpz_t elem;
  mpz_init(elem);
  
  for(int i = 1; i <= rows; i++){
    gmp_matrix_get_elem(elem, i, chainNumber, _Hbasis[dim]);
    elemli = mpz_get_si(elem);
    elemi = elemli;
    coeffVector.push_back(elemi);
    //printf("coeff: %d \n", coeffVector.at(i-1));
  }
  
  mpz_clear(elem);
  return coeffVector;  
}

gmp_matrix* ChainComplex::getBasis(int dim, int basis)
{
  if(dim > -2 && dim < 5 && basis == 2) return _codH[dim+1];
  if(dim < 0 || dim > 4) return NULL;
  if(basis == 0) return create_gmp_matrix_identity(getBasisSize(dim, 0));
  else if(basis == 1) return _kerH[dim];
  else if(basis == 3) return _Hbasis[dim];
  else return NULL;
}

void ChainComplex::getBasisChain(std::map<Cell*, int, Less_Cell>& chain, 
				 int num, int dim, int basis)
{
  gmp_matrix* basisMatrix;
  if(basis == 0) basisMatrix = getBasis(dim, 0);
  else if(basis == 1) basisMatrix = getBasis(dim, 1);
  else if(basis == 2) basisMatrix = getBasis(dim, 2);
  else if(basis == 3) basisMatrix = getBasis(dim, 3);
  else return;

  chain.clear();
  if(dim < 0 || dim > 4) return;
  if(basisMatrix == NULL
     || (int)gmp_matrix_cols(basisMatrix) < num) return;

  int rows = gmp_matrix_rows(basisMatrix);

  int elemi;
  long int elemli;
  mpz_t elem;
  mpz_init(elem);

  for(citer cit = firstCell(dim); cit != lastCell(dim); cit++){
    Cell* cell = cit->first;
    int index = cit->second;
    gmp_matrix_get_elem(elem, index, num, basisMatrix); 
    elemli = mpz_get_si(elem);
    elemi = elemli;
    if(elemli != 0){
      std::map<Cell*, int, Less_Cell > subCells;
      cell->getCells(subCells);
      for(Cell::citer it = subCells.begin(); it != subCells.end(); it++){
	Cell* subCell = (*it).first;
	int coeff = (*it).second;
	chain[subCell] = coeff*elemi; 
      }
    }
  }
  mpz_clear(elem);
}

int ChainComplex::getBasisSize(int dim, int basis)
{
    gmp_matrix* basisMatrix;
    if(basis == 0 && _HMatrix[dim] != NULL){
      return gmp_matrix_cols(_HMatrix[dim]);
    }
    else if(basis == 1) basisMatrix = getBasis(dim, 1);
    else if(basis == 2) basisMatrix = getBasis(dim, 2);
    else if(basis == 3) basisMatrix = getBasis(dim, 3);
    else return 0;
    
    if(basisMatrix != NULL) return gmp_matrix_cols(basisMatrix);
    else return 0; 
}


int ChainComplex::getTorsion(int dim, int num)
{
  if(dim < 0 || dim > 4) return 0;
  if(_Hbasis[dim] == NULL 
     || (int)gmp_matrix_cols(_Hbasis[dim]) < num) return 0;
  if(_torsion[dim].empty() 
     || (int)_torsion[dim].size() < num) return 1;
  else return _torsion[dim].at(num-1);
}


HomologySequence::HomologySequence(ChainComplex* subcomplex, 
				   ChainComplex* complex, 
				   ChainComplex* relcomplex)
{
  _subcomplex = subcomplex;
  _complex = complex;
  _relcomplex = relcomplex;

  mpz_t elem;
  mpz_init_set_si(elem, -1);  

  for(int i = 0; i < 4; i++){
    
    //printf("Dimension %d. \n", i);
    _Ic_sub[i] = NULL;
    _Ic_rel[i] = NULL;

    _Dh[i] = NULL;
    _invDh[i] = NULL;

    _Jh[i] = NULL;
    _Ih[i] = NULL;
    _invJh[i] = NULL;
    _invIh[i] = NULL;

    mpz_t one;
    mpz_init_set_si(one, 1);
    if(_complex->getBasisSize(i, 0) > 0 
       && _subcomplex->getBasisSize(i, 0) > 0){
      _Ic_sub[i] = create_gmp_matrix_zero(_complex->getBasisSize(i, 0), 
					  _subcomplex->getBasisSize(i, 0));
      //printf("rows %d, cols %d. \n", _complex->getBasisSize(i, 0),
      //	     _subcomplex->getBasisSize(i, 0));
      for(ChainComplex::citer cit = _complex->firstCell(i);
	  cit != _complex->lastCell(i); cit++){
	Cell* cell = cit->first;
	int row = cit->second;
	int col = _subcomplex->cellIndex(cell);
	//printf("row %d, col %d. \n", row, col);
	if(col != 0) gmp_matrix_set_elem(one, row, col, _Ic_sub[i]);
      }
    }

    if(_complex->getBasisSize(i, 0) > 0 
       && _relcomplex->getBasisSize(i, 0) > 0){
      _Ic_rel[i] = create_gmp_matrix_zero(_complex->getBasisSize(i, 0), 
					  _relcomplex->getBasisSize(i, 0));
      //printf("rows %d, cols %d. \n", _complex->getBasisSize(i, 0), 
      for(ChainComplex::citer cit = _complex->firstCell(i);
	  cit != _complex->lastCell(i); cit++){
	Cell* cell = cit->first;
	int row = cit->second;
	int col = _relcomplex->cellIndex(cell);
	//printf("row %d, col %d. \n", row, col);
	if(col != 0) gmp_matrix_set_elem(one, row, col, _Ic_rel[i]);
      }
    }
    mpz_clear(one);

    if(_Ic_sub[i] != NULL 
       && _complex->getBasisSize(i, 3) > 0 
       && _subcomplex->getBasisSize(i, 3) > 0){
      gmp_matrix* IH = copy_gmp_matrix(_Ic_sub[i], 1, 1, 
				       gmp_matrix_rows(_Ic_sub[i]), 
				       gmp_matrix_cols(_Ic_sub[i]));
      gmp_matrix_right_mult(IH, _subcomplex->getBasis(i, 3));
      _Ih[i] = createIncMap(IH, _complex->getBasis(i, 3)); 
    }
    if(_Ic_sub[i] != NULL 
       && _complex->getBasisSize(i, 3) > 0 
       && _subcomplex->getBasisSize(i, 3) > 0){
      gmp_matrix* IH = copy_gmp_matrix(_Ic_sub[i], 1, 1, 
				       gmp_matrix_rows(_Ic_sub[i]), 
				       gmp_matrix_cols(_Ic_sub[i]));
      gmp_matrix_transp(IH);
      gmp_matrix_right_mult(IH, _complex->getBasis(i, 3));
      _invIh[i] = createIncMap(IH, _subcomplex->getBasis(i, 3)); 
    }

    if(_Ic_rel[i] != NULL 
       && _complex->getBasisSize(i, 3) > 0 
       && _relcomplex->getBasisSize(i, 3) > 0){
      gmp_matrix* JH = copy_gmp_matrix(_Ic_rel[i], 1, 1, 
				       gmp_matrix_rows(_Ic_rel[i]), 
				       gmp_matrix_cols(_Ic_rel[i]));
      gmp_matrix_transp(JH);
      gmp_matrix_right_mult(JH, _complex->getBasis(i, 3));
      _Jh[i] = createIncMap(JH, _relcomplex->getBasis(i, 3)); 
    }
    if(_Ic_rel[i] != NULL 
       && _complex->getBasisSize(i, 3) > 0 
       && _relcomplex->getBasisSize(i, 3) > 0){
      gmp_matrix* JH = copy_gmp_matrix(_Ic_rel[i], 1, 1, 
				       gmp_matrix_rows(_Ic_rel[i]), 
				       gmp_matrix_cols(_Ic_rel[i]));
      gmp_matrix_right_mult(JH, _relcomplex->getBasis(i, 3));
      _invJh[i] = createIncMap(JH, _complex->getBasis(i, 3)); 
    }

    //printMatrix(_Ih[i]);
    //printMatrix(_invIh[i]);
    //printMatrix(_Jh[i]);
    //printMatrix(_invJh[i]);
    
    if(i > 0 && _relcomplex->getBasisSize(i, 3) > 0 
       && _subcomplex->getBasisSize(i-1, 3) > 0
       && _complex->getBoundaryOp(i) != NULL){
      gmp_matrix* JDIH = copy_gmp_matrix(_Ic_sub[i-1], 1, 1, 
					 gmp_matrix_rows(_Ic_sub[i-1]), 
					 gmp_matrix_cols(_Ic_sub[i-1]));
      gmp_matrix_transp(JDIH);
      gmp_matrix_right_mult(JDIH, _complex->getBoundaryOp(i));
      gmp_matrix_right_mult(JDIH, _Ic_rel[i]);
      gmp_matrix_right_mult(JDIH, _relcomplex->getBasis(i, 3));
      _Dh[i] = createIncMap(JDIH, _subcomplex->getBasis(i-1, 3));
    }

    if(i > 0 && _relcomplex->getBasisSize(i, 3) > 0 
       && _subcomplex->getBasisSize(i-1, 3) > 0
       && _complex->getBoundaryOp(i) != NULL){
      gmp_matrix* JDIH = copy_gmp_matrix(_Ic_rel[i], 1, 1, 
					 gmp_matrix_rows(_Ic_rel[i]), 
					 gmp_matrix_cols(_Ic_rel[i]));
      gmp_matrix_transp(JDIH);
      gmp_matrix* bd = _complex->getBoundaryOp(i);
      gmp_matrix_transp(bd);
      gmp_matrix_right_mult(JDIH, bd);
      gmp_matrix_transp(bd);
      gmp_matrix_right_mult(JDIH, _Ic_sub[i-1]);
      gmp_matrix_right_mult(JDIH, _subcomplex->getBasis(i-1, 3));
      _invDh[i] = createIncMap(JDIH, _relcomplex->getBasis(i, 3));
    }
    
    //printMatrix(_Dh[i]);
    //printMatrix(_invDh[i]);

  }
  for(int i = 0; i < 3; i++){
    blockHBasis(_Dh[i+1], _invIh[i], _subcomplex, i);
    blockHBasis(_Ih[i], _invJh[i], _complex, i);
    blockHBasis(_Jh[i], _invDh[i], _relcomplex, i);
  }
  
}

HomologySequence::~HomologySequence() 
{
  for(int i = 0; i < 4; i++){
    destroy_gmp_matrix(_Ic_sub[i]);
    destroy_gmp_matrix(_Ic_rel[i]);
    destroy_gmp_matrix(_Ih[i]);
    destroy_gmp_matrix(_Jh[i]);
    destroy_gmp_matrix(_invIh[i]);
    destroy_gmp_matrix(_invJh[i]);
    destroy_gmp_matrix(_Dh[i]);
    destroy_gmp_matrix(_invDh[i]);
  }
}

//i: a->b  : aBasis = bBasis*incMap
gmp_matrix* HomologySequence::createIncMap(gmp_matrix* domBasis, 
					   gmp_matrix* codBasis) 
{
  if(domBasis == NULL || codBasis == NULL){
    printf("ERROR: null matrix given. \n");
    return NULL;
  }
  
  int rows = gmp_matrix_rows(domBasis);
  int cols = gmp_matrix_cols(domBasis);
  if(rows < cols || rows == 0 || cols == 0) return NULL;
 
  rows = gmp_matrix_rows(codBasis);
  cols = gmp_matrix_cols(codBasis);
  if(rows < cols || rows == 0 || cols == 0) return NULL;

  gmp_matrix* temp = codBasis;  
  codBasis = copy_gmp_matrix(temp, 1, 1,
			     gmp_matrix_rows(temp),
			     gmp_matrix_cols(temp));
  // inv(U)*A*inv(V) = S
  gmp_normal_form* normalForm 
    = create_gmp_Smith_normal_form(codBasis, INVERTED, INVERTED);
  
  mpz_t elem;
  mpz_init(elem);

  for(int i = 1; i <= cols; i++){
    gmp_matrix_get_elem(elem, i, i, normalForm->canonical);
    if(mpz_cmp_si(elem,0) == 0){
      destroy_gmp_normal_form(normalForm);
      return NULL;
    }
  }
  
  gmp_matrix_left_mult(normalForm->left, domBasis); 
  gmp_matrix* LB = copy_gmp_matrix(domBasis, 1, 1, 
				   gmp_matrix_cols(codBasis), 
				   gmp_matrix_cols(domBasis));
  destroy_gmp_matrix(domBasis);

  rows = gmp_matrix_rows(LB);
  cols = gmp_matrix_cols(LB);

  mpz_t divisor;
  mpz_init(divisor);
  mpz_t remainder;
  mpz_init(remainder);
  mpz_t result;
  mpz_init(result);

  for(int i = 1; i <= rows; i++){
    gmp_matrix_get_elem(divisor, i, i, normalForm->canonical);
    for(int j = 1; j <= cols; j++){
      gmp_matrix_get_elem(elem, i, j, LB);
      mpz_cdiv_qr(result, remainder, elem, divisor);
      if(mpz_cmp_si(remainder, 0) == 0){
        gmp_matrix_set_elem(result, i, j, LB);
      }
      else return NULL;
    }
  }

  gmp_matrix_left_mult(normalForm->right, LB);

  mpz_clear(elem);
  mpz_clear(divisor);
  mpz_clear(result);
  destroy_gmp_normal_form(normalForm);
  return LB;
}


gmp_matrix* HomologySequence::removeZeroCols(gmp_matrix* matrix) // FIXME
{
  mpz_t elem;
  mpz_init(elem);

  int rows = gmp_matrix_rows(matrix);
  int cols = gmp_matrix_cols(matrix);

  std::vector<int> zcols;

  for(int j = 1; j <= cols; j++){
    bool zcol = true;
    for(int i = 1; i <= rows; i++){
      gmp_matrix_get_elem(elem, i, j, matrix);
      if(mpz_cmp_si(elem, 0) != 0){
	zcol = false;
	break;
      }
    }
    if(zcol) zcols.push_back(j);
  }
  if(zcols.empty()) return matrix;
  
  gmp_matrix* newMatrix = create_gmp_matrix_zero(rows, cols-zcols.size());
  if(cols-zcols.size() == 0) return newMatrix;

  int k = 0;
  for(int j = 1; j <= cols; j++){
    if(zcols.size()-1 < k) break;
    if(zcols.at(k) == j) { k++; continue; }
    for(int i = 1; i <= rows; i++){
      gmp_matrix_get_elem(elem, i, j, matrix);
      gmp_matrix_set_elem(elem, i, j-k, newMatrix);
    }
  }

  destroy_gmp_matrix(matrix);
  mpz_clear(elem);
  return newMatrix;
}

void HomologySequence::blockHBasis(gmp_matrix* block1T, 
				   gmp_matrix* block2T, 
				   ChainComplex* complex, int dim)
{
  if(block1T == NULL && block2T == NULL) return;

  gmp_matrix* Hbasis = complex->getBasis(dim, 3);
  
  if(block1T == NULL && block2T != NULL){
    gmp_matrix_right_mult(Hbasis, block2T);
    return;
  }
  if(block1T != NULL && block2T == NULL){
    gmp_matrix_right_mult(Hbasis, block1T);
    return;
  }

  int rows = gmp_matrix_rows(Hbasis);
  int cols = gmp_matrix_cols(Hbasis);
  gmp_matrix* temp1 = copy_gmp_matrix(Hbasis, 1, 1, rows, cols); 
  gmp_matrix* temp2 = copy_gmp_matrix(Hbasis, 1, 1, rows, cols);

  gmp_matrix_right_mult(temp1, block1T); 
  gmp_matrix_right_mult(temp2, block2T);
  
  
  temp1 = removeZeroCols(temp1);
  temp2 = removeZeroCols(temp2);
  //printMatrix(temp1);
  //printMatrix(temp2);

  int bcol = gmp_matrix_cols(temp1);
  
  mpz_t elem;
  mpz_init(elem);
  for(int i = 1; i <= rows; i++){
    for(int j = 1; j <= cols; j++){
      if(j <= bcol) gmp_matrix_get_elem(elem, i, j, temp1);
      else gmp_matrix_get_elem(elem, i, j-bcol, temp2);
      gmp_matrix_set_elem(elem, i, j, Hbasis);
    }
  }

  //printMatrix(Hbasis);
  mpz_clear(elem);
  destroy_gmp_matrix(temp1);
  destroy_gmp_matrix(temp2);
}


#endif
