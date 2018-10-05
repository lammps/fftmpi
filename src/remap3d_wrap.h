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

/* C interface to fftMPI library, 3d Remap functions */

#include <mpi.h>
#include <stdint.h>

#include "fftdata.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

/* C interface */

void remap3d_create(MPI_Comm, void **);
void remap3d_create_fortran(MPI_Comm, void **);
void remap3d_destroy(void *);
void remap3d_set(void *, char *, int);
void remap3d_setup(void *,
                   int, int, int, int, int, int, int, int, int, int, int, int,
                   int, int, int, int *, int *);
void remap3d_setup_fortran(void *,
                           int, int, int, int, int, int, 
                           int, int, int, int, int, int,
                           int, int, int, int *, int *);
void remap3d_remap(void *, FFT_SCALAR *, FFT_SCALAR *, 
                   FFT_SCALAR *, FFT_SCALAR *);

#ifdef __cplusplus
}
#endif
