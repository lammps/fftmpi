/* C interface to fftMPI library, 2d Remap functions */

#include <mpi.h>
#include <stdint.h>

#include "fftdata.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

/* C interface */

void remap2d_create(MPI_Comm, void **);
void remap2d_create_fortran(MPI_Comm, void **);
void remap2d_destroy(void *);
void remap2d_set(void *, char *, int);
void remap2d_setup(void *,
                   int, int, int, int, int, int, int, int,
                   int, int, int, int *, int *);
void remap2d_setup_fortran(void *,
                           int, int, int, int, int, int, int, int,
                           int, int, int, int *, int *);
void remap2d_remap(void *, FFT_SCALAR *, FFT_SCALAR *, 
                   FFT_SCALAR *, FFT_SCALAR *);

#ifdef __cplusplus
}
#endif
