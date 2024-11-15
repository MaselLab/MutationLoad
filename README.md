Whole-genome forwards-time simulation code for the Masel Lab at the University of Arizona, which is able
to simulate the whole genome efficiently via the use of linkage blocks, with recombination only occurring
at hotspots. Specific information about mutations is then recovered with tskit to limit both time and space
complexity.

Version dependencies:
  tskit -- 0.5.8
  gsl -- 2.8
