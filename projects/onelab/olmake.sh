#!/bin/bash

export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab
export MATHEX_DIR=$GMSH_DIR/contrib/MathEx

echo ""
echo "ONELAB: compile a console loader"
g++ -g -I $ONELAB_DIR -I $COMMON_DIR -I $MATHEX_DIR $MATHEX_DIR/mathex.cpp $COMMON_DIR/StringUtils.cpp $COMMON_DIR/onelabUtils.cpp $ONELAB_DIR/myOS.cpp $ONELAB_DIR/OnelabMessage.cpp $ONELAB_DIR/OnelabParser.cpp   $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/loader.cpp -o loader


