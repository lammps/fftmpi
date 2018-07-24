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
SRC2D =	fft2d.cpp remap2d.cpp fft2d_wrap.cpp memory.cpp error.cpp
INC2D =	ffttype.h fft2d.h remap2d.h pack2d.h fft2d_wrap.h memory.h error.h
OBJ2D = $(SRC2D:.cpp=.o)

LIB3D =   libfft3dmpi.a
SHLIB3D = libfft3dmpi.so
SRC3D =	fft3d.cpp remap3d.cpp fft3d_wrap.cpp memory.cpp error.cpp
INC3D =	ffttype.h fft3d.h remap3d.h pack3d.h fft3d_wrap.h memory.h error.h
OBJ3D = $(SRC3D:.cpp=.o)

# memory, 1d FFT, precision settings

MEM_INC =	-DFFT_MEMALIGN=64

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

EXTRA_INC =     $(MEM_INC) $(FFT_INC) $(PRECISION)

# targets

lib:	lib2d lib3d

shlib:	shlib2d shlib3d

all:	shlib lib

help:
	@echo 'make                default = lib, KISS 1d fft, double precision'
	@echo 'make lib            build static fftMPI libs: 2d and 3d'
	@echo 'make shlib          build shared fftMPI libs: 2d and 3d'
	@echo 'make all            build 4 fftMPI libs: lib and shlib'
	@echo 'make lib2d          build static 2d fftMPI lib
	@echo 'make lib3d          build static 3d fftMPI lib
	@echo 'make shlib2d        build shared 2d fftMPI lib
	@echo 'make shlib3d        build shared 3d fftMPI lib
	@echo 'make ... fft=fftw   build with FFTW3 (FFTW), FFTW2, MKL 1d ffts'
	@echo 'make ... p=single   build for single precision'
	@echo 'make clean          remove all *.o files'
	@echo 'make clean-all      remove *.o and lib files'
	@echo 'make tar            create a tarball, 2 levels up'

lib2d:
	$(MAKE) clean
	$(MAKE) static2d fft=$(fft) p=$(p)

lib3d:
	$(MAKE) clean
	$(MAKE) static3d fft=$(fft) p=$(p)

shlib2d:
	$(MAKE) clean
	$(MAKE) shared2d fft=$(fft) p=$(p)

shlib3d:
	$(MAKE) clean
	$(MAKE) shared3d fft=$(fft) p=$(p)

static2d:	$(OBJ2D)
		$(ARCHIVE) $(ARCHFLAGS) $(LIB2D) $(OBJ2D)

static3d:	$(OBJ3D)
		$(ARCHIVE) $(ARCHFLAGS) $(LIB3D) $(OBJ3D)

shared2d:	$(OBJ2D)
		$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) -o $(SHLIB2D) $(OBJ2D)

shared3d:	$(OBJ3D)
		$(CC) $(CCFLAGS) $(SHFLAGS) $(SHLIBFLAGS) -o $(SHLIB3D) $(OBJ3D)

clean:
	@rm -f *.o *.pyc

clean-all:
	@rm -f *.o *.pyc lib*.a lib*.so

tar:
	cd ../..; tar cvf cslib.tar cslib/README cslib/LICENSE \
		cslib/doc cslib/src cslib/test

# rules

%.o:%.cpp
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $<
