#include <cstdio>
#include <vector>

void generate2dEdgeClosureFull(std::vector<std::vector<int> > &closure,
			       std::vector<int> &closureRef,
			       int order, int nNod, bool serendip)
{
  closure.clear();
  closure.resize(2*nNod);
  closureRef.resize(2*nNod);
  int shift = 0;
  for (int corder = order; corder>=0; corder -= (nNod == 3 ? 3 : 2)) {
    if (corder == 0) {
      for (int r = 0; r < nNod ; r++){
        closure[r].push_back(shift);
        closure[r+nNod].push_back(shift);
      }
      break;
    }
    for (int r = 0; r < nNod ; r++){
      for (int j = 0; j < nNod; j++) {
        closure[r].push_back(shift + (r + j) % nNod);
        closure[r + nNod].push_back(shift + (r - j + 1 + nNod) % nNod);
      }
    }
    shift += nNod;
    int n = nNod*(corder-1);
    for (int r = 0; r < nNod ; r++){
      for (int j = 0; j < n; j++){
        closure[r].push_back(shift + (j + (corder - 1) * r) % n);
        closure[r + nNod].push_back(shift + (n - j - 1 + (corder - 1) * (r + 1)) % n);
      }
    }
    shift += n;
    if (serendip) break;
  }
  for (int r = 0; r < nNod*2 ; r++) {
    //closure[r].type = polynomialBasis::getTag(TYPE_LIN, order);
    closureRef[r] = 0;
  }
}

void generate2dEdgeClosure(std::vector<std::vector<int> > &closure, int order, int nNod)
{
  closure.clear();
  closure.resize(2*nNod);
  for (int j = 0; j < nNod ; j++){
    closure[j].push_back(j);
    closure[j].push_back((j+1)%nNod);
    closure[nNod+j].push_back((j+1)%nNod);
    closure[nNod+j].push_back(j);
    for (int i=0; i < order-1; i++){
      closure[j].push_back( nNod + (order-1)*j + i );
      closure[nNod+j].push_back(nNod + (order-1)*(j+1) -i -1);
    }
    //closure[j].type = closure[nNod+j].type = polynomialBasis::getTag(TYPE_LIN, order);
  }
}


int main(void){
  std::vector<std::vector<int> > closure;
  std::vector<std::vector<int> > closureFull;
std::vector<int> closureRef;

  int order = 3;
  int nNode = 3;

  generate2dEdgeClosureFull(closureFull, closureRef, order, nNode, false);

  printf("ClusureFull:\n");

  for(int i = 0; i < closureFull.size(); i++){
    for(int j = 0; j < closureFull[i].size(); j++)
      printf("%d\t", (closureFull[i])[j]);

    printf("\n");
  }

  printf("\n");

  printf("ClosureRef:\n");
  
  for(int i = 0; i < closureRef.size(); i++)
    printf("%d\t", closureRef[i]);
  printf("\n");


  generate2dEdgeClosure(closure, order, nNode);

  printf("Clusure:\n");

  for(int i = 0; i < closure.size(); i++){
    for(int j = 0; j < closure[i].size(); j++)
      printf("%d\t", (closure[i])[j]);

    printf("\n");
  }

  printf("\n");
}
