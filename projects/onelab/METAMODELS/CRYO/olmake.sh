export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile metamodel 'cryosurgery'" 
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp metamodel.cpp -o cryosurgery
