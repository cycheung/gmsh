#!/bin/bash

for i in {1..5}
do
    for j in {1..5}
    do
        ./swavev mesh.msh 5 $i
    done
done

for i in {1..5}
do
    for j in {1..5}
    do
        ./swavev mesh.msh 5 $i slow
    done
done