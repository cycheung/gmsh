#!/bin/bash

for i in {1..5}
do
    ./swavev mesh.msh $i $i
done

for i in {1..5}
do
    ./swavev mesh.msh $i $i slow
done