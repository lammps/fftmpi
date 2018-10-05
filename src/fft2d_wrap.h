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

/* C interface to fftMPI library, 3d FFT functions */

#include <mpi.h>
#include <stdint.h>

#include "fftdata.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

/* C interface */

void fft2d_create(MPI_Comm, int, void **);
void fft2d_create_fortran(MPI_Fint, int, void **);
void fft2d_destroy(void *);

void fft2d_set(void *, const char *, int);
int fft2d_get_int(void *, const char *);
double fft2d_get_double(void *, const char *);
int64_t fft2d_get_int64(void *, const char *);
char *fft2d_get_string(void *, const char *, int *);
  int *fft2d_get_int_vector(void *, const char *, int *);
double *fft2d_get_double_vector(void *, const char *, int *);

void fft2d_setup(void *, int, int,
                 int, int, int, int, int, int, int, int,
                 int, int *, int *, int *);
void fft2d_setup_fortran(void *, int, int,
                         int, int, int, int, int, int, int, int,
                         int, int *, int *, int *);
void fft2d_setup_memory(void *, FFT_SCALAR *, FFT_SCALAR *);

void fft2d_compute(void *, FFT_SCALAR *, FFT_SCALAR *, int);
void fft2d_only_1d_ffts(void *, FFT_SCALAR *, int);
void fft2d_only_remaps(void *, FFT_SCALAR *, FFT_SCALAR *, int);
void fft2d_only_one_remap(void *, FFT_SCALAR *, FFT_SCALAR *, int, int);

void fft2d_tune(void *, int, int,
                int, int, int, int, int, int, int, int,
                int, int *, int *, int *,
                int, int, double, int);
void fft2d_tune_fortran(void *, int, int,
                        int, int, int, int, int, int, int, int,
                        int, int *, int *, int *,
                        int, int, double, int);

#ifdef __cplusplus
}
#endif
