#include "Gmsh.h"
#include "GmshMessage.h"
#include "Context.h"
#include "GModel.h"
#include "LuaBindings.h"

#include "DgC0PlateSolver.h"
#include "PView.h"
#include "PViewData.h"
#include "groupOfElements.h"
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "IPState.h"
int main (int argc, char* argv[])
{

  if (argc != 2){
    printf("usage : dgplate input_file_name (.dat or .lua) \n");
    return -1;
  }
  GmshInitialize(argc, argv);
  GmshSetOption("General","Terminal",1.);
  // get file extension
  std::string fileName = argv[1];
  size_t ext_pos = fileName.find_last_of('.');
  std::string ext(fileName,ext_pos+1,fileName.size());
  std::string dat("dat");
  std::string lua("lua");

  if(ext==dat){
    // globals are still present in Gmsh

    // instanciate a solver
    DgC0PlateSolver mySolver (1000);

    // read some input file
    mySolver.readInputFile(argv[1]);
    //mySolver.createInterfaceElement();
    // create the interface element after reading the mesh file put this somewhere else ??

    // solve the problem
    switch(mySolver.getScheme()){
      case DgC0PlateSolver::StaticLinear :
        Msg::Info("Static linear resolution");
        mySolver.solve();
        break;
      case DgC0PlateSolver::StaticNonLinear :
        Msg::Info("Static non linear resolution");
        mySolver.solveSNL();
        break;
      default : Msg::Error("The scheme of resolution seems to be missing");
    }
  }
  else if(ext==lua){
    #if defined(HAVE_LUA)
    binding *b = binding::instance();
    partDomain::registerBindings(b);
    dgPartDomain::registerBindings(b);
    dgLinearShellDomain::registerBindings(b);
    materialLaw::registerBindings(b);
    linearElasticLawPlaneStress::registerBindings(b);
    linearElasticLawPlaneStressWithFracture::registerBindings(b);
    DgC0PlateSolver::registerBindings(b);
    GmshMergeFile(argv[1]);
    #else
    Msg::Error("Lua is not correctly installed on your computer. Please install it properly before used .lua file or use .dat file")
    #endif
  }
  else Msg::Error("Bad extension for input file. Please use a file with an extension .dat or .lua");
  // stop gmsh
  GmshFinalize();
}
