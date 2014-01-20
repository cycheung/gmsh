#ifndef _NODESOLUTION_H_
#define _NODESOLUTION_H_

#include <string>
#include <complex>
#include <map>

#include "Mesh.h"
#include "PViewDataGModel.h"

/**
   @class NodeSolution
   @brief TODO
 */

class NodeSolution{
 private:
  PViewDataGModel* pView;

 public:
   NodeSolution(void);
  ~NodeSolution(void);

  void addNodeValue(size_t step,
                    double time,
                    const Mesh& mesh,
                    std::map<const MVertex*, std::complex<double> >& data);

  void write(std::string fileName) const;
};

#endif
