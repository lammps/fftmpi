/* C and Fortran interfaces to FFT library */
/* see library.cpp and library_fortran.cpp for wrapper methods */

#include <mpi.h>
#include <stdint.h>

#include "fftdata.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

/* C interface */

void fft2d_create(MPI_Comm, int, void **);
void fft2d_create_fortran(MPI_Comm, int, void **);
void fft2d_create_no_mpi(int, void **);
void fft2d_destroy(void *);
void fft2d_set(void *, const char *, int);
void *fft2d_get(void *, const char *);
void fft2d_setup(void *, int, int,
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

void remap2d_create(MPI_Comm, void **);
void remap2d_create_fortran(MPI_Comm, void **);
void remap2d_create_no_mpi(void **);
void remap2d_destroy(void *);
void remap2d_set(void *, char *, int);
void remap2d_setup(void *,
                   int, int, int, int, int, int, int, int,
                   int, int, int, int *, int *);
void remap2d_remap(void *, FFT_SCALAR *, FFT_SCALAR *, 
                   FFT_SCALAR *, FFT_SCALAR *);

#ifdef __cplusplus
}
#endif
