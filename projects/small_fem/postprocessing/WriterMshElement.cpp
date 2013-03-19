#include "WriterMsh.h"

using namespace std;

void WriterMsh::writeHeader(void) const{
  *out << "$MeshFormat" << endl
       << "2.2 0 8" << endl
       << "$EndMeshFormat" << endl;
}

void WriterMsh::writeNodes(void) const{
  *out << "$Nodes" << endl
       << N << endl;

  for(unsigned int i = 0; i < N; i++){
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

  for(unsigned int i = 0; i < E; i++){
    *out << (*element)[i]->getNum()        << " "
	 << (*element)[i]->getTypeForMSH()
	 << " 2 1 1" << " ";
           // 2 Tags --> (1 physical entity, 1 elementary geometry)

    const unsigned int M = (*element)[i]->getNumVertices();
    MElement* myElement = const_cast<MElement*>((*element)[i]);

    for(unsigned int j = 0; j < M; j++)
      *out << myElement->getVertex(j)->getNum() << " ";

    *out << endl;
  }

  *out << "$EndElements" << endl;
}
