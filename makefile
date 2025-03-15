# Makefile for fgt
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create
# the file make.inc, which overrides the defaults below (which are
# for ubunutu linux/gcc system).

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
CLINK = -lstdc++
FLINK = $(CLINK)

FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -w
# -pg -no-pie is for profiling
#FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -fcx-limited-range -pg -no-pie

# put this in your make.inc if you have FFTW>=3.3.5 and want thread-safe use...
#CXXFLAGS += -DFFTW_PLAN_SAFE
# FFTW base name, and math linking...
FFTWNAME = fftw3
# linux default is fftw3_omp, since 10% faster than fftw3_threads...
#FFTWOMPSUFFIX = omp
FFTWOMPSUFFIX = threads
LIBS := -lm $(CLINK)

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp
#OMP = OFF

LBLAS = -lblas -llapack

FGT_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FGT_INSTALL_DIR = ${HOME}/lib
endif

ifeq ($(PREFIX_FINUFFT),)
	FINUFFT_INSTALL_DIR = ${HOME}/lib
	FINUFFT_INCLUDE_DIR = ${HOME}/include
endif

# absolute path of this makefile, ie FGT's top-level directory...
FGT = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# For your OS, override the above by placing make variables in make.inc
-include make.inc

INCL = -I$(FINUFFT_INCLUDE_DIR)
# here /usr/include needed for fftw3.f "fortran header"... (JiriK: no longer)
FFLAGS := $(FFLAGS) -I/usr/include $(INCL)

DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

# single-thread total list of math and FFTW libs (now both precisions)...
# (Note: finufft tests use LIBSFFT; spread & util tests only need LIBS)
LIBSFFT := -l$(FFTWNAME) -l$(FFTWNAME)f $(LIBS)

# multi-threaded libs & flags, and req'd flags (OO for new interface)...
ifneq ($(OMP),OFF)
  CXXFLAGS += $(OMPFLAGS)
  CFLAGS += $(OMPFLAGS)
  FFLAGS += $(OMPFLAGS)
  MFLAGS += $(MOMPFLAGS) -DR2008OO
  OFLAGS += $(OOMPFLAGS) -DR2008OO
  LIBS += $(OMPLIBS)
  ifneq ($(MINGW),ON)
    ifneq ($(MSYS),ON)
# omp override for total list of math and FFTW libs (now both precisions)...
      #LIBSFFT := -l$(FFTWNAME) -l$(FFTWNAME)_$(FFTWOMPSUFFIX) -l$(FFTWNAME)f -l$(FFTWNAME)f_$(FFTWOMPSUFFIX) $(LIBS)
    endif
  endif
endif

LIBNAME=$(PREFIX_LIBNAME)
ifeq ($(LIBNAME),)
	LIBNAME=libfgt
endif
ifeq ($(MINGW),ON)
  DYNLIB = lib/$(LIBNAME).dll
else
  DYNLIB = lib/$(LIBNAME).so
endif

DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)
# absolute path to the .so, useful for linking so executables portable...
ABSDYNLIB = $(FGT)$(DYNLIB)

FINUFFTLIBNAME = libfinufft
LFINUFFTLINKLIB = -lfinufft

LLINKLIB = $(subst lib, -l, $(LIBNAME))

FINUFFTSTATICLIB = $(FINUFFT_INSTALL_DIR)/$(FINUFFTLIBNAME).a

#
# Note: the static library is used for DYLIBS, so that fmm3d
# does not get bundled in with the fmm3dbie dynamic library
#
LIBS := $(LIBSFFT)
DYLIBS := -L$(FINUFFT_INSTALL_DIR) $(LFINUFFTLINKLIB) $(LIBSFFT)
F2PYDYLIBS += -L$(FINUFFT_INSTALL_DIR) $(LFINUFFTLINKLIB) $(LIBSFFT)

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)

#
# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/prini_new.o \
	$(COM)/hkrand.o \
	$(COM)/dlaran.o \
	$(COM)/cumsum.o \
	$(COM)/fmmcommon2d.o \
	$(COM)/pts_tree_nd.o \
	$(COM)/tree_routs_nd.o \
	$(COM)/besseljs3d.o \
	$(COM)/legeexps.o \
	$(COM)/chebexps.o \
	$(COM)/polytens.o \
	$(COM)/voltab2d.o \
	$(COM)/voltab3d.o \
	$(COM)/tree_data_routs_nd.o \
	$(COM)/tensor_prod_routs_nd.o \
	$(COM)/lapack_f77.o \
	$(COM)/tree_vol_coeffs_nd.o \
	$(COM)/fgtterms.o

# point Gauss transform objects
PFGT = src/pfgt
PFGTOBJS = $(PFGT)/pfgt.o \
	$(PFGT)/pfgt_direct.o \
	$(PFGT)/pfgt_nufftrouts.o

# box Gauss transform objects
BFGT = src/bfgt
BFGTOBJS = $(BFGT)/boxfgt.o \
	$(BFGT)/bfgt_volrouts.o \
	$(BFGT)/bfgt_pwrouts.o \
	$(BFGT)/bfgt_local.o


# Test objects
OBJS = $(COMOBJS) $(PFGTOBJS) $(BFGTOBJS)



.PHONY: usage lib install test-static test-dyn python

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for fgt. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test-static - compile and run validation tests"
	@echo "  make test-dyn - test successful installation by validation tests linked to dynamic library"
	@echo "  make python - compile and test python interfaces using python"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "-------------------------------------------------------------------------"

#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

#
# build the library...
#
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif

$(STATICLIB): $(OBJS)
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/

$(DYNAMICLIB): $(OBJS)
	$(FC) -shared -Wl,-rpath,$(FINUFFT_INSTALL_DIR) $(FFLAGS) $(DYLIBS) $(OBJS) -o $(DYNAMICLIB)
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/

install: $(STATICLIB) $(DYNAMICLIB)
	echo $(FGT_INSTALL_DIR)
	mkdir -p $(FGT_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FGT_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FGT_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FGT_INSTALL_DIR)/
	@echo "Make sure to include " $(FGT_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FGT_INSTALL_DIR)  " "$(LLINKLIB) " -L"$(FINUFFT_INSTALL_DIR)  " "$(LFINUFFTLINKLIB)


#
# testing routines
#
test-static: $(STATICLIB)  test/pfgt-static test/bfgt-static
	cd test/pfgt; ./int2-pfgt
	cd test/bfgt; ./int2-bfgt

test-dyn: $(DYNAMICLIB)  test/pfgt-dyn test/bfgt-dyn
	cd lib; ../test/pfgt/int2-pfgt
	cd lib; ../test/bfgt/int2-bfgt

test/pfgt-static:
	$(FC) $(FFLAGS) test/pfgt/test_pfgt_all.f -o test/pfgt/int2-pfgt lib-static/$(STATICLIB) $(FINUFFTSTATICLIB) $(LIBS)


test/bfgt-static:
	$(FC) $(FFLAGS) test/bfgt/test_boxfgt_all.f -o test/bfgt/int2-bfgt lib-static/$(STATICLIB) $(FINUFFTSTATICLIB) $(LIBS)


#
# Linking test files to dynamic libraries
#

test/pfgt-dyn:
	$(FC) $(FFLAGS) test/pfgt/test_pfgt_all.f $(ABSDYNLIB) -o test/pfgt/int2-pfgt


test/bfgt-dyn:
	$(FC) $(FFLAGS) test/bfgt/test_boxfgt_all.f $(ABSDYNLIB) -o test/bfgt/int2-bfgt

#
# build the python bindings/interface
#
python: $(STATICLIB)
	cd python && export FGTND_LIBS='$(LIBS)' && pip install -e .

#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/pfgt/int2-pfgt
	rm -f test/bfgt/int2-bfgt

objclean:
	rm -f $(OBJS) $(TOBJS)
	rm -f test/pfgt/*.o
	rm -f test/bfgt/*.o