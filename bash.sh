#!/bin/bash

# run Makefile
make
# set library path for shared libraries
LD_LIBRARY_PATH=/home/wmawass/gsl/lib
export LD_LIBRARY_PATH
# run mutationload program with arguments
./mutationload 10000 5000 2 50 23 0.001 0.001 1 1 0 24 20000 1 0.5
./mutationload 10000 5000 2 50 23 0.001 0.005 1 1 0 24 20000 1 0.5
./mutationload 10000 5000 2 50 23 0.001 0.010 1 1 0 24 20000 1 0.5
./mutationload 10000 5000 2 50 23 0.001 0.015 1 1 0 24 20000 1 0.5
./mutationload 10000 5000 2 50 23 0.001 0.020 1 1 0 24 20000 1 0.5
./mutationload 10000 5000 2 50 23 0.001 0.001 1 1 0 24 20000 1 0.2
./mutationload 10000 5000 2 50 23 0.001 0.005 1 1 0 24 20000 1 0.2
./mutationload 10000 5000 2 50 23 0.001 0.010 1 1 0 24 20000 1 0.2
./mutationload 10000 5000 2 50 23 0.001 0.015 1 1 0 24 20000 1 0.2
./mutationload 10000 5000 2 50 23 0.001 0.020 1 1 0 24 20000 1 0.2
