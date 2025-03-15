# makefile overrides
# OS:       macOS
# Compiler: gfortran 14.X
# OpenMP:   enabled
# BLAS:     framework acceralate
#
# NOTE for user:
#           Check gfortran version number 
#

CC=gcc-14
CXX=g++-14
FC=gfortran-14
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -w

ifeq ($(PREFIX_FINUFFT),)
    FINUFFT_INSTALL_DIR=${HOME}/lib
    FINUFFT_INSTALL_DIR_INC=${HOME}/include
endif

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

LBLAS=-framework accelerate

LIBS += -L/usr/local/lib -L/opt/homebrew/lib
