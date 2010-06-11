//
// C++ Interface: terms
//
// Description: Base class to manage a field on element (All field are managed by element in full Dg so displacement is also
//                                                       considered like a field on element
//
// The base class contains archiving data (which is not dependent of the field)
// Author:  <Gauthier BECKER>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "elementField.h"

void elementField::updateFileName(){
  // add a numero to file name
  std::ostringstream oss;
  oss << ++numfile;
  std::string s = oss.str();
  // cut filename and its extension
  size_t ext_pos;
  if(numfile==1)
    ext_pos = fileName.find_last_of('.');
  else{
    oss.clear();
    oss << numfile-1;
    std::string snum = oss.str();
    ext_pos = fileName.find_last_of(snum);
  }
  std::string newname = fileName;
  newname.resize(ext_pos);
  ext_pos = fileName.find_last_of('.');
  std::string ext(fileName,ext_pos+1,fileName.size());
  fileName  = newname+s+"."+ext;
  // create initialization of file
  FILE *fp = fopen(fileName.c_str(), "w");
  fprintf(fp, "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n");
  fclose(fp);
}

void elementField::createFile(){
  if(view){
    // creation of file to store displacement (delete the previous one if exist)
    size_t ext_pos = fileName.find_last_of('.');
    std::string newname = fileName;
    newname.resize(ext_pos);
    std::string rfn = "rm "+newname+"*";
    system(rfn.c_str());
    FILE *fp = fopen(fileName.c_str(), "w");
    fprintf(fp, "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n");
    fclose(fp);
  }
}

elementField::elementField(const std::string &fnn, const uint32_t fms, const int ncomp, const dataType dt, const std::string fin,
                           const bool view_=true) : numfile(0), fileName(fnn), fmaxsize(fms), view(view_),
                                                      totelem(0), numcomp(ncomp), type(dt), fieldname(fin){
  this->createFile();
  switch(type){
    case ElementData :
      dataname = "ElementData";
      break;
    case ElementNodeData :
      dataname = "ElementNodeData";
      break;
  }
}

void elementField::buildView(const std::vector<DGelasticField> &elasticFields,const double time,
                              const int nstep, const bool binary=false){
  if(view){
    // test file size (and create a new one if needed)
    uint32_t fileSize;
    FILE *fp = fopen(fileName.c_str(), "r");
    if(!fp){
      Msg::Error("Unable to open file '%s'", fileName.c_str());
      return;
    }
    fseek (fp, 0, SEEK_END);
    fileSize = (uint32_t) (ftell(fp));
    if(fileSize > fmaxsize) this->updateFileName();
    fclose(fp);
    fp = fopen(fileName.c_str(), "a");

    // compute the number of element
    if(binary) Msg::Warning("Binary write not implemented yet");
    fprintf(fp, "$%s\n",dataname.c_str());
    fprintf(fp, "1\n\"%s\"\n",fieldname.c_str());
    fprintf(fp, "1\n%.16g\n", time);
    fprintf(fp, "3\n%d\n%d\n%d\n", nstep, numcomp, totelem);
    std::vector<double> fieldData;
    if(type == ElementNodeData){
      for (unsigned int i = 0; i < elasticFields.size(); ++i)
        for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it){
          MElement *ele = *it;
          int numv = ele->getNumVertices();
          fieldData.resize(numcomp*numv);
          this->get(ele,fieldData);
          fprintf(fp, "%d %d",ele->getNum(),numv);
          for(int i=0;i<numv;i++)
            for(int j=0;j<numcomp;j++)
              fprintf(fp, " %.16g",fieldData[i+j*numv]);
          fprintf(fp,"\n");
        }
    }
    else if(type == ElementData){
      for (unsigned int i = 0; i < elasticFields.size(); ++i)
        for (groupOfElements::elementContainer::const_iterator it = elasticFields[i].g->begin(); it != elasticFields[i].g->end(); ++it){
          MElement *ele = *it;
          fieldData.resize(numcomp);
          this->get(ele,fieldData);
          fprintf(fp, "%d",ele->getNum());
          for(int j=0;j<numcomp;j++)
            fprintf(fp, " %.16g",fieldData[j]);
          fprintf(fp,"\n");
        }
    }
    fprintf(fp, "$End%s\n",dataname.c_str());
    fclose(fp);
  }
  else Msg::Warning("No displacement view created because the variable view is set to false for displacement field\n");
}
