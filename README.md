This repository contains code to do entanglement clustering for two 
models. The Manuscript for this project is at

https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.023212

Note that the machine learning algorithm used in this code is
UMAP (https://github.com/lmcinnes/umap). The UMAP API is the same
as other manifold learning algorithms in the SKLearn library.
For details of the algorithm see the appendix in the MS.

See below for a description of the dynamic programming algorithm used for
the toric code simulations.

REQUIREMENTS:
  Python 3, numpy, matplotlib, umap (https://github.com/lmcinnes/umap), 
  sklearn, and Julia 0.6

FREE FERMION SIMULATIONS:
  The file simulate.py runs a simulation of the free fermion model described
  in the MS to calculate snap shots of the swap operator.
  The file vmcsim.py contains the backend simulation utilities, including
  the construction of the Slater determinant wave functions and the 
  Monte Carlo step functions. Note that I use the SMW trick for calculating
  determinat ratios. 

  The UMAP machine learning step takes place in unsupervised.py.
  error_analysis.py contains the k-means clustering of the UMAP output
  for accuracy calculations.

TORIC CODE SIMULATIONS:
  simulate_replica.jl runs a simulation of the 2D toric code as described
  in the MS. To improve the efficiency of these simulations, I use a 
  dynamic programming approach to the simulation. I precompute an array
  the same shape as the lattice counting the number of spinons up and to the
  left of each site. 
