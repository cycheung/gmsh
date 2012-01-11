Installation de ONELAB

>> cd XXX/gmsh/projects/onelab
>> olbuild.sh

Le répertoire SOLVERS contient les interfaces vers les solveurs
getdp, elmer...

Ces interfaces sont indépendantes de tout (méta)modèle.
Elles s'attendent à recevoir lors de l'"analyze" des fichiers de données
(getdp: name.pro, elmer:name.sif_onelab)
les informations nécessaires à l'exécutaion de la simulation:

Gmsh/MshFileName
Solver/1ModelName
Solver/1ResolutionChoices (pour getdp)
Solver/2PostOperationChoices (pour getdp)
Solver/9Output files
Solver/MshFileName  (qui est en général différent de Gmsh/MshFileName pour elmer)

L'utilisation des clients solveurs génériques (getdp, elmer)
est limité à la chaîne de résolution standard:
"maillage/préprocessing/calcul/postprocessing"

L'objectif du projet onelab est aussi de permettre l'élaboration de modèles
multicodes plus élaborés, appelés métamodèles.
Pour plus de détails, lire 

./METAMODELS/CRYO/README.txt
./METAMODELS/CORE/README.txt

