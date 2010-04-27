//
// Description : BAR elasticity solver, element space function enriched on tagged vertex
//
//
// Author:  <Eric Bechet>::<Boris Sedji>,  01/2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef _BARELASTICITY_SOLVER_H_
#define _BARELASTICITY_SOLVER_H_

#include "elasticitySolver.h"
#include "simpleFunction.h"
#include "BeamParam.h"
#include "BeamLin.h"

class BarElasticitySolver : public elasticitySolver
{
  protected :
    // map containing the tag of vertex and enriched status
    std::set<int > _TagEnrichedVertex;
    // enriched comp
    std::set<int>  _EnrichComp;
    // simple multiplying function enrichment to enrich the space function
    simpleFunction<double> *_funcEnrichment;
    // beam definitions
    std::vector <BeamLin> _Beams ;

  public :

    BarElasticitySolver(int tag) : elasticitySolver(tag) {}
    ~BarElasticitySolver();
    // create a GModel and determine de dimension of mesh in meshFileName
    virtual void setMesh(const std::string &meshFileName);
    // system solve, read the .dat file, fill tagEnrichedVertex, fill funcEnrichment, solve
    virtual void solve();
    virtual PView *buildDisplacementView(const std::string &postFileName);

    void BeamsIntegration();
};

#endif
