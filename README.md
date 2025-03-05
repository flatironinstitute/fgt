# FGT
Fast Gauss transform (FGT) for discrete and continuous sources in two and three dimensions

# Introduction

This repository contains a simple implementation of the new version of the fast Gauss transform that evaluates the discrete and continuous Gauss transforms:

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


