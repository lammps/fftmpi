/* ----------------------------------------------------------------------
   fftMPI - library for computing 3d/2d FFTs in parallel
   http://fftmpi.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level fftMPI directory.
------------------------------------------------------------------------- */

#ifndef FFT_FFTTYPE_H
#define FFT_FFTTYPE_H

// FFT_PRECISION = 1 for single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 for double-precision complex (8-byte real, 8-byte imag)

#include "fftdata.h"

#ifdef FFT_SINGLE
#define FFT_PRECISION 1
#else
#define FFT_PRECISION 2
#endif

// set default fftw library to FFT_FFTW3

#ifdef FFT_FFTW
#define FFT_FFTW3
#endif

// -------------------------------------------------------------------------

// data types for single-precision complex

#if FFT_PRECISION == 1

#if defined(FFT_MKL)
#include "mkl_dfti.h"
typedef float _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_SINGLE

#elif defined(FFT_FFTW2)
#if defined(FFTW_SIZE)
#include "sfftw.h"
#else
#include "fftw.h"
#endif
typedef FFTW_COMPLEX FFT_DATA;

#elif defined(FFT_FFTW3)
#include "fftw3.h"
typedef fftwf_complex FFT_DATA;
#define FFTW_API(function)  fftwf_ ## function

#else

// use a stripped down version of KISS as default fft

#ifndef FFT_KISS
#define FFT_KISS
#endif

#define kiss_fft_scalar float
typedef struct {
  kiss_fft_scalar re;
  kiss_fft_scalar im;
} FFT_DATA;

struct kiss_fft_state;
typedef struct kiss_fft_state* kiss_fft_cfg;
#endif

// -------------------------------------------------------------------------

// data types for double-precision complex

#elif FFT_PRECISION == 2

#if defined(FFT_MKL)
#include "mkl_dfti.h"
typedef double _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_DOUBLE

#elif defined(FFT_FFTW2)
#if defined(FFTW_SIZE)
#include "dfftw.h"
#else
#include "fftw.h"
#endif
typedef FFTW_COMPLEX FFT_DATA;

#elif defined(FFT_FFTW3)
#include "fftw3.h"
typedef fftw_complex FFT_DATA;
#define FFTW_API(function)  fftw_ ## function

#else

// use a stripped down version of KISS as default fft

#ifndef FFT_KISS
#define FFT_KISS
#endif

#define kiss_fft_scalar double
typedef struct {
  kiss_fft_scalar re;
  kiss_fft_scalar im;
} FFT_DATA;

struct kiss_fft_state;
typedef struct kiss_fft_state* kiss_fft_cfg;
#endif

// -------------------------------------------------------------------------

#else
#error "FFT_PRECISION needs to be either 1 (=single) or 2 (=double)"
#endif

#endif
