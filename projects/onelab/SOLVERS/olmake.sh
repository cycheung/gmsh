export COMMON_DIR=$GMSH_DIR/Common
export ONELAB_DIR=$GMSH_DIR/projects/onelab

echo ""
echo "ONELAB: compile onelab interface client for 'elmer'" 
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp ./elmer.cpp $ONELAB_DIR/main.cpp -o elmer

echo "ONELAB: compile onelab interface client for 'elast'" 
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp ./elast.cpp $ONELAB_DIR/main.cpp -o elast

exit 

echo "ONELAB: compile onelab interface client for 'getdp'" 
g++ -I $ONELAB_DIR -I $COMMON_DIR $COMMON_DIR/OS.cpp $ONELAB_DIR/OnelabClients.cpp $ONELAB_DIR/OnelabMessage.cpp ./getdp.cpp $ONELAB_DIR/main.cpp -o getdp

