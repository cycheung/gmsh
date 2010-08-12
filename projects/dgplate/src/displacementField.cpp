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

#include "displacementField.h"

// constructor
displacementField::displacementField(dofManager<double> *pas, std::vector<DGelasticField> &elas,
                                     const int nc, const int field, const std::vector<Dof> &archiving,
                                     const bool view_, const std::string filen
                                                                         ) : pAssembler(pas), fullDg(elas[0].getFormulation()),
                                                                            _field(field),
                                                                            elementField(filen,1000000,nc,elementField::ElementNodeData,view_)
{
  pAssembler->getFixedDof(fixedDof);
  long int totelem=0;
  // key depends on formulation (change this !!) //TODO ??
  if(!fullDg){ // cG/dG
    // loop on element to initialize ufield
    for(int i=0;i<elas.size();i++){
      for(groupOfElements::vertexContainer::const_iterator it = elas[i].g->vbegin(); it != elas[i].g->vend(); ++it){
        MVertex *ver = *it;
        std::vector<double> u;
        for(int j=0;j<numcomp;j++) u.push_back(0.);
        umap.insert(std::pair<long int,std::vector<double> >(ver->getNum(),u));
      }
    }
    for (unsigned int i = 0; i < elas.size(); ++i)
      for (groupOfElements::elementContainer::const_iterator it = elas[i].g->begin(); it != elas[i].g->end(); ++it)
        totelem++;
   // copy node to archive
   for(int i=0;i<archiving.size();i++)
     varch.insert(std::pair<Dof,long int>(archiving[i],archiving[i].getEntity()));
  }
  else{ // full Dg
    for(int i=0;i<elas.size();i++){
      for(groupOfElements::elementContainer::const_iterator it = elas[i].g->begin(); it != elas[i].g->end(); ++it){
        totelem++;
        MElement *ele = *it;
        int nbval = numcomp*ele->getNumVertices();
        std::vector<double> u;
        for(int j=0;j<nbval;j++) u.push_back(0.);
        umap.insert(std::pair<long int,std::vector<double> >(ele->getNum(),u));
      }
    }
    // Archiving find the element linked to the node
    for(int i=0;i<archiving.size();i++){
      long int nodenum = archiving[i].getEntity();
      int comp,field,num;
      DgC0PlateDof::getThreeIntsFromType(archiving[i].getType(),comp,field,num);
      bool flagout = false;
      // find the element
      for(int k=0;k<elas.size();k++)
        for(groupOfElements::elementContainer::const_iterator it = elas[k].g->begin(); it != elas[k].g->end(); ++it){
          MElement *ele = *it;
          for(int j=0;j<ele->getNumVertices();j++){
            if(ele->getVertex(j)->getNum() == nodenum){
              varch.insert(std::pair<Dof,long int >(Dof(ele->getNum(),DgC0PlateDof::createTypeWithThreeInts(comp,field,j)),nodenum));
              flagout = true;
              break;
            }
          }
          if(flagout) break;
        }
    }
  }
  this->setTotElem(totelem); // TODO FIX IT HOW ??
  this->buildView(elas,0.,0,"displacement",-1,false);
}

void displacementField::update(){
  if(!fullDg){ // Cg/Dg formulation
    for(std::map<long int, std::vector<double> >::iterator it=umap.begin(); it!=umap.end();++it){
      long int ent = (*it).first;
      std::vector<double> u;
      //u.resize(numcomp);
      double du;
      for(int j=0;j<numcomp;j++){
        Dof R=Dof(ent,DgC0PlateDof::createTypeWithThreeInts(j,_field));
        pAssembler->getDofValue(R,du);
        bool in=false;
        for(std::vector<Dof>::iterator itD=fixedDof.begin(); itD<fixedDof.end();++itD) if(*itD==R){in=true;break;}
        if(in) u.push_back(du);
        else   u.push_back((*it).second[j]+du);
      }
      umap[ent]=u;
      u.clear();
    }
  }
  else{ // fullDg formulation
    for(std::map<long int, std::vector<double> >::iterator it=umap.begin(); it!=umap.end();++it){
      long int ent = (*it).first;
      std::vector<double> u;
      int usize=(*it).second.size();
      //u.resize(usize);
      double du;
      int numnodes = usize/numcomp;
      for(int i=0;i<numnodes;i++){
        for(int j=0;j<numcomp;j++){
          Dof R=Dof(ent,DgC0PlateDof::createTypeWithThreeInts(j,_field,i));
          bool in =false;
          pAssembler->getDofValue(R,du);
          for(std::vector<Dof>::iterator itD=fixedDof.begin(); itD<fixedDof.end();++itD) if(*itD==R){in=true;break;}
          if(in)
            u.push_back(du);
          else
            u.push_back((*it).second[i*numcomp+j]+du);
        }
      }
      umap[ent]=u;
      u.clear();
    }
  }
}

void displacementField::get(Dof &D,double &udof){
  long int ent=D.getEntity();
  int field,comp,num;
  DgC0PlateDof::getThreeIntsFromType(D.getType(),comp,field,num);
  if(field!=_field){Msg::Warning("try to get displacement for a non displacement Dof\n"); udof=0;}
  else udof = umap.find(ent)->second[numcomp*num+comp];
}

void displacementField::get(MVertex *ver,std::vector<double> &udofs){
  long int ent;
  if(!fullDg){
    ent = ver->getNum();
    udofs = umap.find(ent)->second;
  }
  else{
    Msg::Error("Impossible to get displacement by vertex for full Dg formulation. Work only with MElement");
  }
}

void displacementField::get(MElement *ele,std::vector<double> &udofs, const int cc){
  long int ent;
  double du;
  if(!fullDg){
    // nb vertices
    int nbvertex = ele->getNumVertices();
    for(int j=0;j<numcomp;j++){
      for(int i=0;i<nbvertex;i++){
        ent = ele->getVertex(i)->getNum();
        Dof D= Dof(ent,DgC0PlateDof::createTypeWithThreeInts(j,_field));
        this->get(D,du);
        udofs.push_back(du);
      }
    }
  }
  else{
    ent = ele->getNum();
    int nbvertex = ele->getNumVertices();
    for(int j=0;j<numcomp;j++){
      for(int i=0;i<nbvertex;i++){
        Dof D= Dof(ent,DgC0PlateDof::createTypeWithThreeInts(j,_field,i));
        this->get(D,du);
        udofs.push_back(du);
      }
    }
  }
}

void displacementField::get(MInterfaceElement* iele, std::vector<double> &udofs){
//  udofs.resize(0);
//  std::vector<double> utemp;
//  int usizem = numcomp*iele->getElem(0)->getNumVertices();
//  utemp.resize(usizem);
  this->get(iele->getElem(0),udofs);
//  for(int i=0;i<usizem;i++)
//    udofs.push_back(utemp[i]);
  if(!(iele->getElem(0)==iele->getElem(1))){// Virtual interface element
    //int usizep=numcomp*iele->getElem(1)->getNumVertices();
    //utemp.resize(usizep);
    this->get(iele->getElem(1),udofs);
    //for(int i=0;i<usizep;i++)
    //  udofs.push_back(utemp[i]);
  }
}

void displacementField::getForPerturbation(MInterfaceElement* iele, const bool minus, Dof &D, double pert, std::vector<double> &udofs){
  double eps = -pert;
  if(!minus){
    this->set(D,eps);
    this->get(iele->getElem(0),udofs);
    this->set(D,pert);
  }
  else{
    this->get(iele->getElem(0),udofs);
    this->set(D,eps);
  }
  if(!(iele->getElem(0)==iele->getElem(1)))// Virtual interface element
    this->get(iele->getElem(1),udofs);
  if(minus)
    this->set(D,pert);
}

void displacementField::updateFixedDof(){
  double u;
  long int ent;
  int type,field,comp,numv;
  for(std::vector<Dof>::iterator itD=fixedDof.begin(); itD<fixedDof.end();++itD){
    pAssembler->getDofValue(*itD,u);
    ent = itD->getEntity();
    type= itD->getType();
    DgC0PlateDof::getThreeIntsFromType(type,comp,field,numv);
    umap.find(ent)->second[numv*numcomp+comp]=u;
  }

}
void displacementField::archiving(const double time){
  FILE *fp;
  double u;
  for(std::map<Dof,long int>::iterator it = varch.begin(); it!=varch.end();++it){
    Dof D = it->first;
    this->get(D,u);
    // write in File
      std::ostringstream oss;
      oss << it->second;
      std::string s = oss.str();
      // component of displacement
      int field,comp,num;
      DgC0PlateDof::getThreeIntsFromType(it->first.getType(),comp,field,num);
      oss.str("");
      oss << comp;
      std::string s2 = oss.str();
      std::string fname = "NodalDisplacement"+s+"comp"+s2+".csv";
      fp = fopen(fname.c_str(),"a");
      fprintf(fp,"%lf;%lf\n",time,u);
      fclose(fp);
  }
}
