# FGT
Fast Gauss transform (FGT) for discrete and continuous sources in two and three dimensions

# Introduction

This repository contains an openmp implementation of the new version of the fast Gauss transform that evaluates the discrete and continuous Gauss transforms:

$$u_i=\sum_{j=1}^N G({\boldsymbol x}_i- {\boldsymbol y}_j;\delta)q_j,
\qquad {\boldsymbol x}_i, {\boldsymbol y}_j \in B, \qquad i=1,\ldots M,$$

and

$$u({\boldsymbol x})=\int_B G({\boldsymbol x}-{\boldsymbol y};\delta)\sigma({\boldsymbol y})d{\boldsymbol y},$$

where $B = [-\frac{1}{2},\frac{1}{2}]^d$ is the 
unit box centered at the origin in $\mathbb{R}^d$. For free space problems,
the Gaussian kernel is given by 

$$G({\boldsymbol x};\delta)=e^{-\frac{\|{\boldsymbol x}\|^2}{\delta}}, $$

while for periodic problems, 

$$G({\boldsymbol x};\delta)= \sum_{{\bf j} \in {\mathbb{Z}^d}} e^{-\frac{\|{\boldsymbol x} + {\bf j}\|^2}{\delta}},$$

where ${\mathbb{Z}}^d$ denotes the $d$-dimensional integer lattice.
For discrete sources, the scheme relies on the nonuniform fast Fourier transform 
(NUFFT) to construct near field plane wave representations. For continuous source 
distributions sampled on adaptive tensor-product grids, we exploit 
the separable structure of the Gaussian kernel to accelerate the 
computation. 

# Installation instructions for point FGT
1. Download the latest version of finufft

2. Compile finufft by "make lib OMP=OFF", i.e., the finufft library
should be compiled in single threaded mode. The current library compiled
with multithreaded finufft may lead to segmentation faults. 

3. Copy libfinufft.a from finufft/lib-static and libfinufft.so from finufft/lib
to the directory that contains your other local library files

4. Copy "finufft.fh" from the directory "~/finufft/include" to the
directory that contains your other local library files

5. Copy over the appropriate ``make.inc.*`` to ``make.inc`` and run ``make
install``

For example: 
* On a mac with gfortran-12: cp make.inc.macos.gnu make.inc && make
install

To verify successful installation of static and dynamic libraries, 
run ``make test-static`` and ``make test-dyn`` respectively. 
The output of the code should have an error of roughly 1e-6 as the last line.

# Callable subroutines

* The driver for discrete summation is src/pfgt/pfgt.f. See the comments in pfgt.f for its
  input and output arguments and test/pfgt/test_pfgt_all.f for examples. 

* The driver for continuous convolution is src/bfgt/boxfgt.f, which requires calling subroutines
vol_tree_mem and vol_tree_build in src/common/tree_vol_coeffs.f first to build the tree. See the
comments in boxfgt.f for its input and output arguments and test/bfgt/test_boxfgt_all.f for examples.

# Citing

If you find FGT useful in your work, please star this repository and cite it and the following. 

```
@article{GJRW2024sirev,
author = {Greengard, Leslie F. and Jiang, Shidong and Rachh, Manas and Wang, Jun},
title = {A New Version of the Adaptive Fast Gauss Transform for Discrete and Continuous Sources},
journal = {SIAM Review},
volume = {66},
number = {2},
pages = {287-315},
year = {2024},
doi = {10.1137/23M1572453},
URL = {https://doi.org/10.1137/23M1572453},
eprint = {https://doi.org/10.1137/23M1572453},
}
```

# Main developers

* Leslie Greengard, Flatiron Institute, Simons Foundation
* Shidong Jiang, Flatiron Institute, Simons Foundation
* Manas Rachh, Flatiron Institute, Simons Foundation
* Jun Wang, Tsinghua University, China
