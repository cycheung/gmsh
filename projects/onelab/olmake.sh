#!/bin/bash

export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile the console (non gui) onelab loader"
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $COMMON_DIR/StringUtils.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp $ONELAB_DIR/OnelabParser.cpp $ONELAB_DIR/SOLVERS/onelab.cpp $ONELAB_DIR/loader.cpp -o loader

echo "ONELAB: compile a simple remote client for testing socket communication" 
g++ -I $ONELAB_DIR -I $COMMON_DIR remote.cpp -o remote
