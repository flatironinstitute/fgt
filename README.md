# FGT — Fast Gauss Transform

Fast Gauss transform (FGT) for discrete and continuous sources in two and three dimensions.

## Introduction

This repository contains an OpenMP implementation of the new version of the fast Gauss transform that evaluates the discrete and continuous Gauss transforms:

$$u_i=\sum_{j=1}^N G({\boldsymbol x}_i- {\boldsymbol y}_j;\delta)q_j,
\qquad {\boldsymbol x}_i, {\boldsymbol y}_j \in B, \qquad i=1,\ldots M,$$

and

$$u({\boldsymbol x})=\int_B G({\boldsymbol x}-{\boldsymbol y};\delta)\sigma({\boldsymbol y})d{\boldsymbol y},$$

where $B = [-\frac{1}{2},\frac{1}{2}]^d$ is the unit box centered at the origin in $\mathbb{R}^d$. For free-space problems, the Gaussian kernel is

$$G({\boldsymbol x};\delta)=e^{-\frac{\|{\boldsymbol x}\|^2}{\delta}},$$

while for periodic problems,

$$G({\boldsymbol x};\delta)= \sum_{{\bf j} \in {\mathbb{Z}^d}} e^{-\frac{\|{\boldsymbol x} + {\bf j}\|^2}{\delta}},$$

where ${\mathbb{Z}}^d$ denotes the $d$-dimensional integer lattice. For discrete sources, the scheme relies on the nonuniform fast Fourier transform (NUFFT) to construct near-field plane-wave representations. For continuous source distributions sampled on adaptive tensor-product grids, we exploit the separable structure of the Gaussian kernel.

## Dependencies

| Component | Required version | Notes |
|---|---|---|
| `gfortran` | ≥ 13 (15 recommended) | Tested with Homebrew `gfortran-15`. |
| `gcc` / `g++` | ≥ 13 | Same toolchain as `gfortran`. |
| FFTW3 | ≥ 3.3.5, single + double precision | `apt install libfftw3-dev`, `brew install fftw`, or `pacman -S mingw-w64-x86_64-fftw`. |
| BLAS / LAPACK | any modern | Apple Accelerate on macOS; OpenBLAS / MKL on Linux. |
| OpenMP runtime | matches compiler | `libgomp` for GCC. |
| Python | ≥ 3.9 | Optional (Python interface). |

`fgt` bundles its FINUFFT dependency as a git submodule pinned to upstream `flatironinstitute/finufft` v2.5.1. The fgt makefile builds it for you — there is no need to download or compile FINUFFT separately.

## Quick install

### macOS (Homebrew, recommended)

```bash
brew install gcc fftw                      # gcc gives you gfortran
git clone --recurse-submodules https://github.com/flatironinstitute/fgt.git
cd fgt
cp make.inc.macos.gnu make.inc             # auto-detects Homebrew gfortran-N
make install OMP=ON                        # builds finufft submodule + fgt
make test-static                           # validates correctness
make perftest                              # validates OMP scaling
```

### Linux (Debian/Ubuntu)

```bash
sudo apt install gfortran libfftw3-dev libopenblas-dev
git clone --recurse-submodules https://github.com/flatironinstitute/fgt.git
cd fgt
cp make.inc.linux.gnu.openblas make.inc
make install OMP=ON
make test-static
```

### Linux (RHEL / Almalinux 8+ / Fedora)

```bash
sudo dnf install gcc-gfortran fftw-devel openblas-devel
git clone --recurse-submodules https://github.com/flatironinstitute/fgt.git
cd fgt
cp make.inc.linux.gnu.openblas make.inc
make install OMP=ON
make test-static
```

### Windows (MSYS2 mingw-w64 shell)

```bash
pacman -S mingw-w64-x86_64-{gcc,gcc-fortran,fftw,openblas,pkgconf,git,make}
git clone --recurse-submodules https://github.com/flatironinstitute/fgt.git
cd fgt
cp make.inc.windows.mingw make.inc
mingw32-make install OMP=ON
mingw32-make test-static
```

If `make install` complains that finufft headers are missing, your clone is missing the submodule. Fix it with:

```bash
git submodule update --init --recursive
```

## The `make.inc` system

`fgt` reads compiler/library settings from `make.inc`, which you create by copying one of the templates:

| Template | Use case |
|---|---|
| `make.inc.macos.gnu` | macOS + Homebrew gfortran + Apple Accelerate. **Recommended on macOS.** Auto-detects the latest Homebrew `gcc-N`/`gfortran-N`. |
| `make.inc.macos.gnu.openblas` | macOS + Homebrew gfortran + OpenBLAS. Note: Homebrew OpenBLAS pulls in Apple `libomp` and conflicts with `libgomp`; the `bfgt` validation test fails with this combo. Prefer Accelerate on Mac. |
| `make.inc.macos.intel` | Intel Mac + Intel oneAPI (`icx`/`ifx`) + MKL. |
| `make.inc.linux.gnu.openblas` | Linux + system GCC + OpenBLAS (`apt install libopenblas-dev libfftw3-dev`, etc.). |
| `make.inc.linux.intel` | Linux + Intel oneAPI + MKL. |
| `make.inc.windows.mingw` | Windows + MSYS2 mingw-w64 (`pacman -S mingw-w64-x86_64-{gcc,gcc-fortran,fftw,openblas}`). |

Override variables (set in `make.inc` or as `make` arguments):

| Variable | Meaning |
|---|---|
| `FC`, `CC`, `CXX` | Compilers. |
| `FFLAGS`, `CFLAGS`, `CXXFLAGS` | Compiler flags. |
| `LBLAS` | BLAS/LAPACK link flags (e.g. `-framework accelerate`, `-lopenblas`, `-mkl`). |
| `FFT_INSTALL_DIR`, `FFT_INCLUDE_DIR` | FFTW install/include dirs. |
| `PREFIX` | Where `make install` copies fgt's libs (default `$(HOME)/lib`). |
| `PREFIX_FINUFFT` | Override path to finufft static lib (defaults to bundled submodule). |
| `PREFIX_FINUFFT_INCLUDE` | Override path to finufft includes (defaults to bundled submodule). |
| `OMP=ON` / `OMP=OFF` | Whether to build fgt with OpenMP (default ON). |

## Build targets

| Target | What it does |
|---|---|
| `make lib` | Builds `lib-static/libfgt.a` and `lib/libfgt.so` (or `.dll` on Windows). |
| `make install` | Copies the libs into `PREFIX`. |
| `make test-static` | Links tests against the static lib and runs them. |
| `make test-dyn` | Same as `test-static` but using the dynamic lib. |
| `make perftest` | Runs an OMP scaling sweep over `OMP_NUM_THREADS = 1, 2, 4, ..., physcores`. Saves CSV to `docs/perf/`. |
| `make python` | Installs the Python wrapper editable into the current environment. |
| `make finufft` | Builds the bundled finufft submodule explicitly (rarely needed; `lib` depends on this). |
| `make clean` / `make objclean` | Removes build outputs. |

## Threading model

`fgt` uses OpenMP at its top level (`pfgt` parallelizes over boxes via `$OMP PARALLEL DO`). Each OMP thread carries its own FINUFFT plan in a per-thread `fftplan(ithd)` array, and FINUFFT itself is built **single-threaded** (the bundled finufft submodule is compiled with `OMP=OFF`). Inside the source we also set `opts%nthreads = 1` as defensive insurance.

**Do not call `pfgt` or `boxfgt` from inside another OpenMP parallel region** — that would oversubscribe and is not tested.

## Performance

Scaling table on Apple M-series, gfortran-15, dim=3, N=10⁶, eps=1e-6, charges only, free-space:

| nthreads | total_s | speedup | efficiency |
|---------:|--------:|--------:|-----------:|
|        1 |   6.567 |    1.00 |       100% |
|        2 |   3.309 |    1.98 |        99% |
|        4 |   1.730 |    3.80 |        95% |
|        8 |   0.957 |    6.86 |        86% |
|       12 |   0.834 |    7.87 |        66% |

The drop in efficiency at the highest core count reflects Apple's mix of performance + efficiency cores (8P + 4E on this CPU). On homogeneous Linux servers expect strong scaling closer to ideal at all power-of-two thread counts up to physical cores.

If your scaling efficiency is much lower, try pinning threads:

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

Raw scaling-sweep CSVs from `make perftest` are saved to `docs/perf/`.

## Python interface

```bash
pip install fgt                         # from PyPI (when available)
# or, from a local clone:
make lib OMP=ON
pip install -e .                        # editable install
```

Minimal example:

```python
import numpy as np
import fgt

rng = np.random.default_rng(0)
ns = 1_000_000
sources = rng.uniform(size=(3, ns))
charges = rng.standard_normal((1, ns))

out = fgt.pfgt(sources, charges=charges, delta=1e-3, eps=1e-6, ifpgh=1)
print(out["pot"].shape)                 # (1, 1000000)
```

The Python wrapper is a thin `ctypes` layer over `libfgt`; the heavy lifting is in Fortran.

## Callable subroutines

* The driver for discrete summation is `src/pfgt/pfgt.f`. See the comments in `pfgt.f` for its input and output arguments and `test/pfgt/test_pfgt_all.f` for examples.
* The driver for continuous convolution is `src/bfgt/boxfgt.f`, which requires calling subroutines `vol_tree_mem` and `vol_tree_build` in `src/common/tree_vol_coeffs.f` first to build the tree. See the comments in `boxfgt.f` for its input and output arguments and `test/bfgt/test_boxfgt_all.f` for examples.

## Citing

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

## Main developers

* Leslie Greengard, Flatiron Institute, Simons Foundation
* Shidong Jiang, Flatiron Institute, Simons Foundation
* Manas Rachh, Flatiron Institute, Simons Foundation
* Jun Wang, Tsinghua University, China
