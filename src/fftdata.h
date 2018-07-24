// FFT data types for single and double precision

#ifndef FFT_FFTDATA_H
#define FFT_FFTDATA_H

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#endif
