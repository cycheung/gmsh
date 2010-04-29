//
// C++ Interface: terms
//
// Description: Insertion of Interface Elements in GModel
//
//
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef GMODELDG_H_
#define GMODELDG_H_
#include "GModel.h"
#include "MInterfaceElement.h"
#include "DgC0PlateSolver.h" // remove if elasticfields is not pass to generateInterfaceElementsOnBoundary
#include "groupOfElements.h" // remove if elasticfields is not pass to generateInterfaceElementsOnBoundary
class GModelWithInterface : public GModel{

  protected :
    // Add Interface element to the GModel I keep the MInterfaceElement because there is no need of a GInterfaceElement
    // I use vector because I don't know how to use set
    std::vector<std::pair<int,MInterfaceElement> > interfaces; //
    std::vector<std::pair<int,MInterfaceElement> > boundinter; // A different vector is used for boundary element of the interface
    std::vector<MInterfaceElement> virtualinterface; // Use to applied boundary condition
    std::vector<std::pair<int,std::vector<MVertex*> > > thetaBound; // Use to store the BC applied on theta (the vertex where are computed the condition are compute after the reading of the mesh )
    // int value needed for fullDg formulation
    //int max_vertex; // max number of vertex of an element in the GModel

    // Store interface elements when I understand this operation
    // void _storeInterfaceElements(std::vector<std::map<int,MInterfaceElement> > & ie );

  public :
    // build function
    GModelWithInterface() : GModel(){};

    // return interface element of a field
    std::vector<MInterfaceElement*> getInterface(const int num_field){
      std::vector<MInterfaceElement*> vie;
      for(int i=0;i<interfaces.size();i++)
        if(interfaces[i].first == num_field) vie.push_back(&interfaces[i].second);
      return vie;
    }
    // return all internal interface element
    std::vector<MInterfaceElement*> getInterface(){
      std::vector<MInterfaceElement*> vie;
      for(int i=0;i<interfaces.size();i++)
        vie.push_back(&interfaces[i].second);
      return vie;
    }

    // Function to generate interface element on boundary
    void generateInterfaceElementsOnBoundary(const int &num_phys, std::vector<DGelasticField> elasticFields);

    // Function create virtual interface element to applied Dirichlet boundary conditions
    void generateVirtualInterfaceElement(const int &num_phys, std::vector<DGelasticField> elasticFields);

    // Return the boundary interfaceElement linked to an elasticField
    std::vector<MInterfaceElement*> getBoundInterface(const int num_field){
      std::vector<MInterfaceElement*> vie;
      for(int i=0;i<boundinter.size();i++)
        if(boundinter[i].first == num_field) vie.push_back(&boundinter[i].second);
      return vie;
    }
    // Return the virtual interfaceElement linked to an elasticField
    std::vector<MInterfaceElement*> getVirtualInterface(){
      std::vector<MInterfaceElement*> vie;
      for(int i=0;i<virtualinterface.size();i++)
        vie.push_back(&virtualinterface[i]);
      return vie;
    }
    void storeInterfaceElement(const int phys,const MInterfaceElement ie){
      interfaces.push_back(std::pair<int,MInterfaceElement>(phys,ie));
    }
    void storeBoundaryInterfaceElement(const int phys,const MInterfaceElement ie){
      boundinter.push_back(std::pair<int,MInterfaceElement>(phys,ie));
    }
    void storeVirtualInterfaceElement(MInterfaceElement ie){
      virtualinterface.push_back(ie);
    }
    void insertThetaBound(const int phys,const std::vector<MVertex*> &v){thetaBound.push_back(std::pair<int,std::vector<MVertex*> >(phys,v));}
    void getThetaBound(const int phys,std::vector<MVertex*> &vv){
      for(int i=0;i<thetaBound.size();i++)
        if(phys==thetaBound[i].first) vv=thetaBound[i].second;
    }
    void getphysBound(std::vector<int> &v){
      v.clear();
      for(int i=0;i<thetaBound.size();i++) v.push_back(thetaBound[i].first);
    }
};


#endif // GMODELDG_H_
