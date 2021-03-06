# Makefile for fftMPI = 2d and 3d parallel FFTs
# type "make help" for options

SHELL = /bin/sh

# ----------------------------------------
# should only need to change this section
# compiler/linker settings
# ----------------------------------------

CC =		mpicxx
CCFLAGS =	-g -O3 
SHFLAGS =	-fPIC
ARCHIVE =       ar
ARCHFLAGS =     -rc
SHLIBFLAGS =	-shared

# files

LIB2D =   libfft2dmpi.a
SHLIB2D = libfft2dmpi.so
SRC2D =	fft2d.cpp remap2d.cpp fft2d_wrap.cpp remap2d_wrap.cpp \
	memory.cpp error.cpp
INC2D =	ffttype.h fft2d.h remap2d.h pack2d.h fft2d_wrap.h remap2d_wrap.h \
	memory.h error.h
OBJ2D = $(SRC2D:.cpp=.o)

LIB3D =   libfft3dmpi.a
SHLIB3D = libfft3dmpi.so
SRC3D =	fft3d.cpp remap3d.cpp fft3d_wrap.cpp remap3d_wrap.cpp \
	memory.cpp error.cpp
INC3D =	ffttype.h fft3d.h remap3d.h pack3d.h fft3d_wrap.h remap2d_wrap.h \
	memory.h error.h
OBJ3D = $(SRC3D:.cpp=.o)

# memory, 1d FFT, precision settings

MEM_INC =	-DFFT_MEMALIGN=64

tbb = 		yes
TBB =		$(shell echo $(tbb) | tr a-z A-Z)

ifeq ($(TBB),NO)
  TBB_INC = -DFFT_INTEL_NO_TBB
endif

fft = 		kiss
FFT =		$(shell echo $(fft) | tr a-z A-Z)
FFT_INC =       -DFFT_$(FFT)

ifeq ($(FFT),MKL)
  FFT_INC = -DFFT_$(FFT) -I$(MKL_INCLUDE)
endif

p = 		DOUBLE
P =       	$(shell echo $(p) | tr a-z A-Z)

ifeq ($(P),SINGLE)
  PRECISION = -DFFT_SINGLE
else
  PRECISION = -DFFT_DOUBLE
endif

EXTRA_INC =     $(MEM_INC) $(TBB_INC) $(FFT_INC) $(PRECISION)

# 1d FFT libraries used for shared fftMPI lib build

ifeq ($(FFT),FFTW2)
  FFT1D = -lfftw
else ifeq ($(FFT),FFTW3)
  ifeq ($p,single)
    FFT1D = -lfftw3f
  else
    FFT1D = -lfftw3
  endif
else ifeq ($(FFT),FFTW)
  ifeq ($p,single)
    FFT1D = -lfftw3f
  else
    FFT1D = -lfftw3
  endif
else ifeq ($(FFT),MKL)
  ifeq ($c,intel)
    FFT1D = -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
  else
    FFT1D = -L$(MKL_LIB) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
  endif
else
  FFT1D =
endif

# targets

shlib:	shlib2d shlib3d

lib:	lib2d lib3d

all:	shlib lib

help:
	@echo 'make                default = shlib, KISS 1d fft, double precision'
	@echo 'make shlib          build 2 shared fftMPI libs: 2d and 3d'
	@echo 'make lib            build 2 static fftMPI libs: 2d and 3d'
	@echo 'make all            build 4 fftMPI libs: lib and shlib'
	@echo 'make shlib2d        build shared 2d fftMPI lib'
	@echo 'make shlib3d        build shared 3d fftMPI lib'
	@echo 'make lib2d          build static 2d fftMPI lib'
	@echo 'make lib3d          build static 3d fftMPI lib'
	@echo 'make ... fft=fftw   build with FFTW3 (FFTW), FFTW2, MKL 1d ffts'
	@echo 'make ... p=single   build for single precision'
	@echo 'make ... tbb=no     build w/out Intel TBB memory allocation'
	@echo 'make clean          remove all *.o files'
	@echo 'make clean-all      remove *.o and lib files'

shlib2d:
	$(MAKE) clean
	$(MAKE) shared2d fft=$(fft) p=$(p)

shlib3d:
	$(MAKE) clean
	$(MAKE) shared3d fft=$(fft) p=$(p)

lib2d:
	$(MAKE) clean
	$(MAKE) static2d fft=$(fft) p=$(p)

lib3d:
	$(MAKE) clean
	$(MAKE) static3d fft=$(fft) p=$(p)

static2d:	$(OBJ2D)
		$(ARCHIVE) $(ARCHFLAGS) $(LIB2D) $(OBJ2D)

static3d:	$(OBJ3D)
		$(ARCHIVE) $(ARCHFLAGS) $(LIB3D) $(OBJ3D)

shared2d:	$(OBJ2D)
		$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) \
		-o $(SHLIB2D) $(OBJ2D) $(FFT1D)

shared3d:	$(OBJ3D)
		$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) \
		-o $(SHLIB3D) $(OBJ3D) $(FFT1D)

clean:
	@rm -f *.o *.pyc

clean-all:
	@rm -f *.o *.pyc lib*.a lib*.so

# rules

%.o:%.cpp
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<
