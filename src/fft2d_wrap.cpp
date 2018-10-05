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

// C interface to fftMPI library, 2d FFT functions

#include <string.h>
#include <stdlib.h>

#include "fft2d_wrap.h"
#include "fft2d.h"

using namespace FFTMPI_NS;

// ----------------------------------------------------------------------
// 2d FFT library calls
// ---------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create an instance of a 2d FFT and return pointer to it
   pass in MPI communicator to run on
------------------------------------------------------------------------- */

void fft2d_create(MPI_Comm comm, int precision, void **ptr)
{
  FFT2d *fft = new FFT2d(comm,precision);
  *ptr = (void *) fft;
}

// ----------------------------------------------------------------------

void fft2d_create_fortran(MPI_Fint fcomm, int precision, void **ptr)
{
  MPI_Comm ccomm = MPI_Comm_f2c(fcomm);
  FFT2d *fft = new FFT2d(ccomm,precision);
  *ptr = (void *) fft;
}

/* ----------------------------------------------------------------------
   destruct an instance of a 2d FFT
------------------------------------------------------------------------- */

void fft2d_destroy(void *ptr)
{
  FFT2d *fft = (FFT2d *) ptr;
  delete fft;
}

/* ----------------------------------------------------------------------
   set an internal flag, before setup() or compute()
------------------------------------------------------------------------- */

void fft2d_set(void *ptr, const char *keyword, int value)
{
  FFT2d *fft = (FFT2d *) ptr;

  if (strcmp(keyword,"collective") == 0) fft->collective = value;
  else if (strcmp(keyword,"exchange") == 0) fft->exchange = value;
  else if (strcmp(keyword,"pack") == 0) fft->packflag = value;
  else if (strcmp(keyword,"memory") == 0) fft->memoryflag = value;
  else if (strcmp(keyword,"scale") == 0) fft->scaled = value;
  else if (strcmp(keyword,"remaponly") == 0) fft->remaponly = value;
}

/* ----------------------------------------------------------------------
   get value of an internal integer variable
------------------------------------------------------------------------- */

int fft2d_get_int(void *ptr, const char *keyword)
{
  FFT2d *fft = (FFT2d *) ptr;

  if (strcmp(keyword,"collective") == 0) return fft->collective;
  else if (strcmp(keyword,"exchange") == 0) return fft->exchange;
  else if (strcmp(keyword,"pack") == 0) return fft->packflag;
  else if (strcmp(keyword,"npfast1") == 0) return fft->npfast1;
  else if (strcmp(keyword,"npfast2") == 0) return fft->npfast2;
  else if (strcmp(keyword,"npslow1") == 0) return fft->npslow1;
  else if (strcmp(keyword,"npslow2") == 0) return fft->npslow2;
  else if (strcmp(keyword,"npbrick1") == 0) return fft->npbrick1;
  else if (strcmp(keyword,"npbrick2") == 0) return fft->npbrick2;
  else if (strcmp(keyword,"ntrial") == 0) return fft->ntrial;
  else if (strcmp(keyword,"npertrial") == 0) return fft->npertrial;

  return -1;
}

/* ----------------------------------------------------------------------
   get value of an internal double variable
------------------------------------------------------------------------- */

double fft2d_get_double(void *ptr, const char *keyword)
{
  FFT2d *fft = (FFT2d *) ptr;

  if (strcmp(keyword,"setuptime") == 0) return fft->setuptime;

  return -1.0;
}

/* ----------------------------------------------------------------------
   get value of an internal int64 variable
------------------------------------------------------------------------- */

int64_t fft2d_get_int64(void *ptr, const char *keyword)
{
  FFT2d *fft = (FFT2d *) ptr;

  if (strcmp(keyword,"memusage") == 0) return fft->memusage;

  return -1;
}

/* ----------------------------------------------------------------------
   get value of an internal string variable
------------------------------------------------------------------------- */

char *fft2d_get_string(void *ptr, const char *keyword, int *len)
{
  FFT2d *fft = (FFT2d *) ptr;

  if (strcmp(keyword,"fft1d") == 0) {
    *len = strlen(fft->fft1d);
    return (char *) fft->fft1d;
  } else if (strcmp(keyword,"precision") == 0) {
    *len = strlen(fft->precision);
    return (char *) fft->precision;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   get pointer to an internal vector of ints variable
------------------------------------------------------------------------- */

int *fft2d_get_int_vector(void *ptr, const char *keyword, int *len)
{
  FFT2d *fft = (FFT2d *) ptr;

  *len = fft->ntrial;

  if (strcmp(keyword,"cflags") == 0) return fft->cflags;
  else if (strcmp(keyword,"eflags") == 0) return fft->eflags;
  else if (strcmp(keyword,"pflags") == 0) return fft->pflags;

  return NULL;
}

/* ----------------------------------------------------------------------
   get pointer to an internal vector of doubles variable
------------------------------------------------------------------------- */

double *fft2d_get_double_vector(void *ptr, const char *keyword, int *len)
{
  FFT2d *fft = (FFT2d *) ptr;

  *len = fft->ntrial;

  if (strcmp(keyword,"tfft") == 0) return fft->tfft;
  else if (strcmp(keyword,"t1d") == 0) return fft->t1d;
  else if (strcmp(keyword,"tremap") == 0) return fft->tremap;
  else if (strcmp(keyword,"tremap1") == 0) return fft->tremap1;
  else if (strcmp(keyword,"tremap2") == 0) return fft->tremap2;
  else if (strcmp(keyword,"tremap3") == 0) return fft->tremap3;

  return NULL;
}

/* ----------------------------------------------------------------------
   create plan for performing a 2d FFT
------------------------------------------------------------------------- */

void fft2d_setup(void *ptr,
                 int nfast, int nslow,
                 int in_ilo, int in_ihi, int in_jlo, int in_jhi,
                 int out_ilo, int out_ihi, int out_jlo, int out_jhi,
                 int permute, int *fftsize_caller,
                 int *sendsize_caller, int *recvsize_caller)
{
  FFT2d *fft = (FFT2d *) ptr;

  int fftsize,sendsize,recvsize;
  fft->setup(nfast,nslow,
             in_ilo,in_ihi,in_jlo,in_jhi,
             out_ilo,out_ihi,out_jlo,out_jhi,
             permute,fftsize,sendsize,recvsize);
  *fftsize_caller = fftsize;
  *sendsize_caller = sendsize;
  *recvsize_caller = recvsize;
}

/* ----------------------------------------------------------------------
   create plan for performing a 2d FFT
   Fortran interface where indices are 1 to N inclusive
------------------------------------------------------------------------- */

void fft2d_setup_fortran(void *ptr,
                         int nfast, int nslow,
                         int in_ilo, int in_ihi, int in_jlo, int in_jhi,
                         int out_ilo, int out_ihi, int out_jlo, int out_jhi,
                         int permute, int *fftsize_caller,
                         int *sendsize_caller, int *recvsize_caller)
{
  FFT2d *fft = (FFT2d *) ptr;

  int fftsize,sendsize,recvsize;
  fft->setup(nfast,nslow,
             in_ilo-1,in_ihi-1,in_jlo-1,in_jhi-1,
             out_ilo-1,out_ihi-1,out_jlo-1,out_jhi-1,
             permute,fftsize,sendsize,recvsize);
  *fftsize_caller = fftsize;
  *sendsize_caller = sendsize;
  *recvsize_caller = recvsize;
}

/* ----------------------------------------------------------------------
   pass in user memory for a 2d remap send/recv
------------------------------------------------------------------------- */

void fft2d_setup_memory(void *ptr, FFT_SCALAR *sendbuf, FFT_SCALAR *recvbuf)
{
  FFT2d *fft = (FFT2d *) ptr;
  fft->setup_memory(sendbuf,recvbuf);
}

/* ----------------------------------------------------------------------
   perform a 2d FFT
------------------------------------------------------------------------- */

void fft2d_compute(void *ptr, FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  FFT2d *fft = (FFT2d *) ptr;
  fft->compute(in,out,flag);
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 2d FFT, no data movement
------------------------------------------------------------------------- */

void fft2d_only_1d_ffts(void *ptr, FFT_SCALAR *in, int flag)
{
  FFT2d *fft = (FFT2d *) ptr;
  fft->only_1d_ffts(in,flag);
}

/* ----------------------------------------------------------------------
   perform all the remaps in a 2d FFT, but no 1d FFTs
------------------------------------------------------------------------- */

void fft2d_only_remaps(void *ptr, FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  FFT2d *fft = (FFT2d *) ptr;
  fft->only_remaps(in,out,flag);
}

/* ----------------------------------------------------------------------
   perform just a single 2d remap operation
------------------------------------------------------------------------- */

void fft2d_only_one_remap(void *ptr, 
                          FFT_SCALAR *in, FFT_SCALAR *out, int flag, int which)
{
  FFT2d *fft = (FFT2d *) ptr;
  fft->only_one_remap(in,out,flag,which);
}

/* ----------------------------------------------------------------------
   tune settings for fastest FFT: collective, exchange, pack flags
------------------------------------------------------------------------- */

void fft2d_tune(void *ptr, 
                int nfast, int nslow,
                int in_ilo, int in_ihi, int in_jlo, int in_jhi,
                int out_ilo, int out_ihi, int out_jlo, int out_jhi, 
                int permute, int *fftsize_caller, 
                int *sendsize_caller, int *recvsize_caller,
                int flag, int niter, double tmax, int tflag)
{
  FFT2d *fft = (FFT2d *) ptr;

  int fftsize,sendsize,recvsize;
  fft->tune(nfast,nslow,
            in_ilo,in_ihi,in_jlo,in_jhi,
            out_ilo,out_ihi,out_jlo,out_jhi,
            permute,fftsize,sendsize,recvsize,
            flag,niter,tmax,tflag);
  *fftsize_caller = fftsize;
  *sendsize_caller = sendsize;
  *recvsize_caller = recvsize;
}

/* ----------------------------------------------------------------------
   tune settings for fastest FFT: collective, exchange, pack flags
   Fortran interface where indices are 1 to N inclusive
------------------------------------------------------------------------- */

void fft2d_tune_fortran(void *ptr, 
                        int nfast, int nslow,
                        int in_ilo, int in_ihi, int in_jlo, int in_jhi,
                        int out_ilo, int out_ihi, int out_jlo, int out_jhi, 
                        int permute, int *fftsize_caller, 
                        int *sendsize_caller, int *recvsize_caller,
                        int flag, int niter, double tmax, int tflag)
{
  FFT2d *fft = (FFT2d *) ptr;

  int fftsize,sendsize,recvsize;
  fft->tune(nfast,nslow,
            in_ilo-1,in_ihi-1,in_jlo-1,in_jhi-1,
            out_ilo-1,out_ihi-1,out_jlo-1,out_jhi-1,
            permute,fftsize,sendsize,recvsize,
            flag,niter,tmax,tflag);
  *fftsize_caller = fftsize;
  *sendsize_caller = sendsize;
  *recvsize_caller = recvsize;
}
