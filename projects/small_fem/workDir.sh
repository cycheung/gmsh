#!/bin/bash

rm -fr build
mkdir build
cd build

cmake -DDEFAULT=0 ..
