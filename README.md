Whole-genome forwards-time simulation code for the Masel Lab at the University of Arizona, which is able
to simulate the whole genome efficiently via the use of linkage blocks, with recombination only occurring
at hotspots. Specific information about mutations is then recovered with tskit to limit both time and space
complexity.

Version dependencies:
  tskit -- 0.5.8
  gsl -- 2.8

To use tskit with the code you can git clone the tskit repo and then checkout to version 0.5.8 

'''
git clone https://github.com/tskit-dev/tskit.git
cd tskit
git checkout 0.5.8
'''

The makefile assumes that you clone tskit on your home directory.