#The variable 'GMSH_DIR' must be set to the home dir of gmsh on your system
#e.g. export GMSH_DIR=$HOME/SOLVERS/gmsh   

export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile the console (non gui) onelab loader"
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp $ONELAB_DIR/loader.cpp -o loader

echo "ONELAB: compile the metamodels"

cd METAMODELS

cd CORE
olmake.sh
cd ..

cd CRYO
olmake.sh
cd ..

cd ..

echo "ONELAB: compile a simple remote client for testing socket communication" 
g++ -I $ONELAB_DIR -I $COMMON_DIR remote.cpp -o remote
