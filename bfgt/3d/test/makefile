
EXEC = int2-fgt

#HOST = osx
HOST=linux-gfortran
#HOST=linux-ifort
#HOST=linux-gfortran-prof
#HOST=linux-gfortran-openmp

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

FINUFFT = /home/shidong/finufft
ABSDYNLIB = $(FINUFFT)/lib/libfinufft.so
# FFTW base name, and math linking...
FFTWNAME = fftw3

LIBSFFT := -l$(FFTWNAME)  $(LIBS)

ifeq ($(HOST),osx)
FC = gfortran
FFLAGS = -O3 -march=native --openmp -funroll-loops -c -w
FLINK = gfortran -w --openmp -o $(EXEC)
FEND = -lopenblas ${LDFLAGS}
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w
FLINK = gfortran -w -o $(EXEC) $(ABSDYNLIB) $(LIBSFFT)
FEND = -lblas -llapack
#FEND = -lopenblas -L/usr/local/opt/openblas/lib 
endif

ifeq ($(HOST),linux-gfortran-prof)
FC = gfortran
FFLAGS = -O3 -march=native -pg -g -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) -pg
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -c --openmp
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -lopenblas -L/usr/local/opt/openblas/lib 
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O1 -g -c -w -xW -qopenmp 
FLINK = ifort -w -qopenmp -o $(EXEC)
WITH_SECOND = 1
endif


SRC = ../src
FGT3D = ../src
UTILS_DIR = ../src/common


.PHONY: all clean list

SOURCES =  test_bfgt3d.f \
  $(SRC)/common/prini_new.f \
  $(SRC)/common/legeexps.f \
  $(SRC)/common/legetens.f \
  $(SRC)/common/voltab3d.f \
  $(SRC)/bfgt3d.f \
  $(SRC)/g3dterms.f \
  $(SRC)/g3drouts.f \
  $(SRC)/fgt2dvolrouts.f \
  $(SRC)/fgt3dvolrouts.f \
  $(SRC)/tree_routs3d.f \
  $(SRC)/common/hkrand.f \
  $(SRC)/common/dlaran.f \
  $(SRC)/tree_vol_coeffs2.f \
#  $(SRC)/tree_vol.f \

ifeq ($(WITH_SECOND),1)
SOURCES += $(SRC)/second-r8.f
endif

OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(EXEC) 2 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



