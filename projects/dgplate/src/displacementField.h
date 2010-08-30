//
// C++ Interface: terms
//
// Description: Class with the displacement field
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef DISPLACEMENTFIELD_H_
#define DISPLACEMENTFIELD_H_
#include "dofManager.h"
#include "DgC0PlateSolver.h"
#include "PView.h"
#include "PViewData.h"
#include <stdint.h>
#include <stdlib.h>
#include "elementField.h"

class displacementField : public elementField{
  protected :
    std::map<long int ,std::vector<double> > umap; // One entry by Dof entity
    dofManager<double> *pAssembler; // To access to component of equation system template this
    std::vector<Dof> fixedDof;
    bool fullDg; // formulation
    int _field;
    std::map<Dof,long int> varch;

  public :
    displacementField(dofManager<double> *pas, std::vector<partDomain*> &elas,
                      const int nc, const int field, const std::vector<Dof> &archiving,
                      const bool =true, const std::string="disp.msh"
                                                    ) ;
    // update all displacement value
    void update();
    // get Operation
    void get(Dof &D,double &udof);
    void get(MVertex *ver,std::vector<double> &udofs); // works only for cG/dG (TODO fullDg case)
    virtual void get(MElement *ele,std::vector<double> &udofs, const int=-1);
    void get(MInterfaceElement* iele, std::vector<double> &udofs);
    void updateFixedDof();
    void set(Dof &D, double &val){
      long int ent = D.getEntity();
      int field,comp,num;
      DgC0PlateDof::getThreeIntsFromType(D.getType(),comp,field,num);
      umap[ent][num*numcomp+comp] += val;
    }
    void getForPerturbation(MInterfaceElement* iele, const bool minus, Dof &D, double pert, std::vector<double> &udofs);
    void archiving(const double time);
};
#endif // DISPLACEMENTFIELD_H_
