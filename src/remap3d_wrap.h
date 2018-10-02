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
