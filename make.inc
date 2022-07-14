# makefile overrides
# OS:       Linux
# Compiler: gfortran 
# OpenMP:   enabled
# BLAS:     openblas
#

LBLAS=-lopenblas
FFLAGS += -DFFTW_PLAN_SAFE