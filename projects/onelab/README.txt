Installation de ONELAB


>> export OL_DIR= $GMSH_DIR/gmsh/projects/onelab
>> cd $OL_DIR
>> olbuild.sh

Le répertoire SOLVERS contient l'interface onelab:

'$OL_DIR/SOLVERS/onelab_client'

1. Utilisation en ligne de commande:

>> cd METAMODELS/CRYO
>> ../../SOLVERS/onelab_client cryo


2. Utilisation avec gmsh comme loader:

Gmsh gui:
>> gmsh cryo.geo
Click  Menu->Tools->Onelab
Enter 'Client Name'=onelab 
Browse for the location $OL_DIR/SOLVERS/getdp
Click on 'onelab_client'

