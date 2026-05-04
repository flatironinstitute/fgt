# makefile overrides
# OS:       macOS (Intel or Apple Silicon)
# Compiler: Homebrew gfortran (>= 13)
# OpenMP:   enabled
# BLAS:     framework Accelerate
# FFTW:     Homebrew (`brew install fftw`)
#
# NOTE for user:
#   * Set FC/CC/CXX to your Homebrew gfortran version, e.g.
#       make FC=gfortran-15 CC=gcc-15 CXX=g++-15
#     The defaults below auto-detect the latest Homebrew gcc on PATH.
#   * fgt requires gfortran >= 13 (we tested 13, 14, 15).

# Auto-detect the Homebrew gfortran/gcc/g++ binary (e.g. gfortran-15).
# The grep filters out gcc-ar-N, gcc-nm-N, gcc-ranlib-N (which match gcc-*
# but are not the compiler).
HOMEBREW_PREFIX := $(shell brew --prefix 2>/dev/null)
GFORTRAN_BIN := $(shell ls $(HOMEBREW_PREFIX)/bin/gfortran-* 2>/dev/null | grep -E '/gfortran-[0-9]+$$' | sort -V | tail -1)
GCC_BIN      := $(shell ls $(HOMEBREW_PREFIX)/bin/gcc-* 2>/dev/null | grep -E '/gcc-[0-9]+$$' | sort -V | tail -1)
GXX_BIN      := $(shell ls $(HOMEBREW_PREFIX)/bin/g++-* 2>/dev/null | grep -E '/g\+\+-[0-9]+$$' | sort -V | tail -1)

ifeq ($(GFORTRAN_BIN),)
$(error Could not find a Homebrew gfortran. Install one: `brew install gcc`. \
        Or set FC/CC/CXX explicitly on the make command line.)
endif

# Use = (not ?=) so we override the makefile's plain `gcc`/`g++`/`gfortran`
# defaults; users can still override on the command line (those take
# precedence over make.inc).
FC  = $(GFORTRAN_BIN)
CC  = $(GCC_BIN)
CXX = $(GXX_BIN)

FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy

# Default install prefix (overridable via `make install PREFIX=...`).
ifeq ($(PREFIX),)
    FGT_INSTALL_DIR=$(HOME)/lib
endif

# By default we link against the bundled finufft submodule, no override needed.

# FFTW from Homebrew.
ifeq ($(PREFIX_FFT),)
    FFTW_ROOT := $(shell brew --prefix fftw)
    FFT_INSTALL_DIR=$(FFTW_ROOT)/lib
endif

ifeq ($(PREFIX_FFT_INCLUDE),)
    FFTW_ROOT := $(shell brew --prefix fftw)
    FFT_INCLUDE_DIR=$(FFTW_ROOT)/include
endif

# OpenMP with Homebrew gcc on macOS.
OMPFLAGS = -fopenmp
OMPLIBS  = -lgomp

# Apple Accelerate for BLAS/LAPACK.
LBLAS=-framework accelerate
