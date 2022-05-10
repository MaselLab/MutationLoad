#!/bin/bash

# run Makefile
make
# set library path for shared libraries
LD_LIBRARY_PATH=/home/wmawass/gsl/lib
export LD_LIBRARY_PATH
# run mutationload program with arguments
./mutationload 100000 10000 2 50 23 0.001 0.010 1 1 0 24 20000 1 0.5 0.8 0.01
#./mutationload 100000 10000 2 50 23 0.001 0.020 1 1 0 24 20000 1 0.5 0.8 0.01
#./mutationload 100000 16000 2 50 23 0.001 0.010 1 1 0 24 20000 1 0.2 0.8 0.01
#./mutationload 100000 16000 2 50 23 0.001 0.020 1 1 0 24 20000 1 0.2 0.8 0.01
#./mutationload 100000 12500 2 50 23 0.001 0.010 1 1 0 24 25000 1 0.5 0.8 0.01
#./mutationload 100000 12500 2 50 23 0.001 0.020 1 1 0 24 25000 1 0.5 0.8 0.01
#./mutationload 100000 20000 2 50 23 0.001 0.010 1 1 0 24 25000 1 0.2 0.8 0.01
#./mutationload 100000 20000 2 50 23 0.001 0.020 1 1 0 24 25000 1 0.2 0.8 0.01
#./mutationload 100000 17500 2 50 23 0.001 0.010 1 1 0 24 35000 1 0.5 0.8 0.01
#./mutationload 100000 17500 2 50 23 0.001 0.020 1 1 0 24 35000 1 0.5 0.8 0.01
