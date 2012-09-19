I. Introduction

Le champ décrit par le fichier "circle.pos" est le domaine de définition de la fonctionnelle étudiée.
On travaille dans un espace normalisé par la magnétisation de saturation $\vec J_S$.
La fonctionnelle à minimiser s'écrit
$$
\Omega(|\vec J|) = J \atanh{J} +1/2 \ln(J^2-1) - H \cdot J + KK |J-J_P|
$$
avec $J=|\vec J|$.
L'énergie magnétique tendant vers l'infini lorsque $|\vec J| \rightarrow 1$,
le domaine de définition est limité à
$$
{ \vec J, |\vec J| < 0.99 }
$$

II. Calcul de la fonctionnelle

La fonctionnelle est évaluée sur "circle.pos" avec le plugin MathEval de Gmsh.

