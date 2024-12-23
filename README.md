# py-exactdiag-project

This is an exact diagonalization library written in C++ with python bindings. I wrote this to better acquaint myself with C++ and implement symmetries in exact diagonalization of quantum systems. 

Given a lattice, model (e.g. spin one-half, real space bosons, real space fermions) and appropriate symmetries, this will generate the state space and matrix representations of operators. It implements an operator algebra with second quantization and generates the appropriate matrix upon request, which one can then diagonalize with numpy and scipy. It can currently handle the following:

-Spin Half
-Real space bosons
-1D translation symmetry
-1D inversion symmetry
-Total Sz conservation (for spin half)

In progress:
-General space group symmetries
-Real space fermions
-Momentum space bosons and fermions
-Lanczos diagonalization methods (so no need to use numpy!)
-Krylov methods for exponentiation

To use this program, see test.py for appropriate use cases. 
