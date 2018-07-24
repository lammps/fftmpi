// C or Fortran style library interface to FFT library

#include <string.h>
#include <stdlib.h>
#include "fft3d_wrap.h"
#include "fft3d.h"
#include "remap3d.h"

using namespace FFTMPI_NS;

// ----------------------------------------------------------------------
// 3d FFT library calls
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create an instance of a 3d FFT and return pointer to it
   pass in MPI communicator to run on
------------------------------------------------------------------------- */

void fft3d_create(MPI_Comm communicator, int precision, void **ptr)
{
  FFT3d *fft = new FFT3d(communicator,precision);
  *ptr = (void *) fft;
}

// ----------------------------------------------------------------------

void fft3d_create_fortran(MPI_Fint fcomm, int precision, void **ptr)
{
  MPI_Comm ccomm = MPI_Comm_f2c(fcomm);
  FFT3d *fft = new FFT3d(ccomm,precision);
  *ptr = (void *) fft;
}

/* ----------------------------------------------------------------------
   create an instance of a 3d FFT and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   initialize MPI if needed
------------------------------------------------------------------------- */

void fft3d_create_no_mpi(int precision, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

  FFT3d *fft = new FFT3d(communicator,precision);
  *ptr = (void *) fft;
}

/* ----------------------------------------------------------------------
   destruct an instance of a 3d FFT
------------------------------------------------------------------------- */

void fft3d_destroy(void *ptr)
{
  FFT3d *fft = (FFT3d *) ptr;
  delete fft;
}

/* ----------------------------------------------------------------------
   set an internal flag, before setup() or compute()
------------------------------------------------------------------------- */

void fft3d_set(void *ptr, const char *keyword, int value)
{
  FFT3d *fft = (FFT3d *) ptr;

  if (strcmp(keyword,"collective") == 0) fft->collective = value;
  else if (strcmp(keyword,"exchange") == 0) fft->exchange = value;
  else if (strcmp(keyword,"pack") == 0) fft->packflag = value;
  else if (strcmp(keyword,"memory") == 0) fft->memoryflag = value;
  else if (strcmp(keyword,"scale") == 0) fft->scaled = value;
  else if (strcmp(keyword,"remaponly") == 0) fft->remaponly = value;
}

/* ----------------------------------------------------------------------
   get value of an internal value, return as pointer to value(s)
   caller must cast the pointer correctly to access the value(s)
------------------------------------------------------------------------- */

void *fft3d_get(void *ptr, const char *keyword)
{
  FFT3d *fft = (FFT3d *) ptr;

  if (strcmp(keyword,"fft1d") == 0) return (void *) fft->fft1d;
  else if (strcmp(keyword,"precision") == 0) return (void *) fft->precision;
  else if (strcmp(keyword,"collective") == 0) return &fft->collective;
  else if (strcmp(keyword,"exchange") == 0) return &fft->exchange;
  else if (strcmp(keyword,"pack") == 0) return &fft->packflag;
  else if (strcmp(keyword,"memusage") == 0) return &fft->memusage;
  else if (strcmp(keyword,"npfast1") == 0) return &fft->npfast1;
  else if (strcmp(keyword,"npfast2") == 0) return &fft->npfast2;
  else if (strcmp(keyword,"npfast3") == 0) return &fft->npfast3;
  else if (strcmp(keyword,"npmid1") == 0) return &fft->npmid1;
  else if (strcmp(keyword,"npmid2") == 0) return &fft->npmid2;
  else if (strcmp(keyword,"npmid3") == 0) return &fft->npmid3;
  else if (strcmp(keyword,"npslow1") == 0) return &fft->npslow1;
  else if (strcmp(keyword,"npslow2") == 0) return &fft->npslow2;
  else if (strcmp(keyword,"npslow3") == 0) return &fft->npslow3;
  else if (strcmp(keyword,"npbrick1") == 0) return &fft->npbrick1;
  else if (strcmp(keyword,"npbrick2") == 0) return &fft->npbrick2;
  else if (strcmp(keyword,"npbrick3") == 0) return &fft->npbrick3;
  else if (strcmp(keyword,"ntrial") == 0) return &fft->ntrial;
  else if (strcmp(keyword,"npertrial") == 0) return &fft->npertrial;
  else if (strcmp(keyword,"setuptime") == 0) return &fft->setuptime;
  else if (strcmp(keyword,"cflags") == 0) return fft->cflags;
  else if (strcmp(keyword,"eflags") == 0) return fft->eflags;
  else if (strcmp(keyword,"pflags") == 0) return fft->pflags;
  else if (strcmp(keyword,"tfft") == 0) return fft->tfft;
  else if (strcmp(keyword,"t1d") == 0) return fft->t1d;
  else if (strcmp(keyword,"tremap") == 0) return fft->tremap;
  else if (strcmp(keyword,"tremap1") == 0) return fft->tremap1;
  else if (strcmp(keyword,"tremap2") == 0) return fft->tremap2;
  else if (strcmp(keyword,"tremap3") == 0) return fft->tremap3;
  else if (strcmp(keyword,"tremap4") == 0) return fft->tremap4;
  else return NULL;
}

/* ----------------------------------------------------------------------
   create plan for performing a 3d FFT
------------------------------------------------------------------------- */

void fft3d_setup(void *ptr,
                 int nfast, int nmid, int nslow,
                 int in_ilo, int in_ihi, int in_jlo, 
                 int in_jhi, int in_klo, int in_khi,
                 int out_ilo, int out_ihi, int out_jlo, 
                 int out_jhi, int out_klo, int out_khi,
                 int permute, int *fftsize_caller, 
                 int *sendsize_caller, int *recvsize_caller)
{
  FFT3d *fft = (FFT3d *) ptr;

  int fftsize,sendsize,recvsize;
  fft->setup(nfast,nmid,nslow,
             in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
             out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
             permute,fftsize,sendsize,recvsize);
  *fftsize_caller = fftsize;
  *sendsize_caller = sendsize;
  *recvsize_caller = recvsize;
}

/* ----------------------------------------------------------------------
   pass in user memory for a 3d remap send/recv
------------------------------------------------------------------------- */

void fft3d_setup_memory(void *ptr, FFT_SCALAR *sendbuf, FFT_SCALAR *recvbuf)
{
  FFT3d *fft = (FFT3d *) ptr;
  fft->setup_memory(sendbuf,recvbuf);
}

/* ----------------------------------------------------------------------
   perform a 3d FFT
------------------------------------------------------------------------- */

void fft3d_compute(void *ptr, FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  FFT3d *fft = (FFT3d *) ptr;
  fft->compute(in,out,flag);
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 3d FFT, no data movement
------------------------------------------------------------------------- */

void fft3d_only_1d_ffts(void *ptr, FFT_SCALAR *in, int flag)
{
  FFT3d *fft = (FFT3d *) ptr;
  fft->only_1d_ffts(in,flag);
}

/* ----------------------------------------------------------------------
   perform all the remaps in a 3d FFT, but no 1d FFTs
------------------------------------------------------------------------- */

void fft3d_only_remaps(void *ptr, FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  FFT3d *fft = (FFT3d *) ptr;
  fft->only_remaps(in,out,flag);
}

/* ----------------------------------------------------------------------
   perform just a single 3d remap operation
------------------------------------------------------------------------- */

void fft3d_only_one_remap(void *ptr,
                          FFT_SCALAR *in, FFT_SCALAR *out, int flag, int which)
{
  FFT3d *fft = (FFT3d *) ptr;
  fft->only_one_remap(in,out,flag,which);
}

/* ----------------------------------------------------------------------
   tune settings for fastest FFT: collective, exchange, pack flags
------------------------------------------------------------------------- */

void fft3d_tune(void *ptr, 
                int nfast, int nmid, int nslow,
                int in_ilo, int in_ihi, int in_jlo, 
                int in_jhi, int in_klo, int in_khi,
                int out_ilo, int out_ihi, int out_jlo, 
                int out_jhi, int out_klo, int out_khi,
                int permute, int *fftsize_caller, 
                int *sendsize_caller, int *recvsize_caller,
                int flag, int niter, double tmax, int tflag)
{
  FFT3d *fft = (FFT3d *) ptr;

  int fftsize,sendsize,recvsize;
  fft->tune(nfast,nmid,nslow,
            in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
            out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
            permute,fftsize,sendsize,recvsize,
            flag,niter,tmax,tflag);
  *fftsize_caller = fftsize;
  *sendsize_caller = sendsize;
  *recvsize_caller = recvsize;
}

// ----------------------------------------------------------------------
// 3d Remap library calls
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create an instance of a 3d Remap and return pointer to it
   pass in MPI communicator to run on
------------------------------------------------------------------------- */

void remap3d_create(MPI_Comm communicator, void **ptr)
{
  Remap3d *remap = new Remap3d(communicator);
  *ptr = (void *) remap;
}

// ----------------------------------------------------------------------

void remap3d_create_fortran(MPI_Fint fcomm, void **ptr)
{
  MPI_Comm ccomm = MPI_Comm_f2c(fcomm);
  Remap3d *remap = new Remap3d(ccomm);
  *ptr = (void *) remap;
}

/* ----------------------------------------------------------------------
   create an instance of a 3d Remap and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   initialize MPI if needed
------------------------------------------------------------------------- */

void remap3d_create_no_mpi(void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

  Remap3d *remap = new Remap3d(communicator);
  *ptr = (void *) remap;
}

/* ----------------------------------------------------------------------
   destruct an instance of a 3d Remap
------------------------------------------------------------------------- */

void remap3d_destroy(void *ptr)
{
  Remap3d *remap = (Remap3d *) ptr;
  delete remap;
}

/* ----------------------------------------------------------------------
   set an internal flag, before setup()
------------------------------------------------------------------------- */

void remap3d_set(void *ptr, char *keyword, int value)
{
  Remap3d *remap = (Remap3d *) ptr;

  if (strcmp(keyword,"collective") == 0) remap->collective = value;
  else if (strcmp(keyword,"pack") == 0) remap->packflag = value;
}

/* ----------------------------------------------------------------------
   create plan for performing a 3d Remap
------------------------------------------------------------------------- */

void remap3d_setup(void *ptr,
                   int in_ilo, int in_ihi, int in_jlo, 
                   int in_jhi, int in_klo, int in_khi,
                   int out_ilo, int out_ihi, int out_jlo, 
                   int out_jhi, int out_klo, int out_khi,
                   int nqty, int permute, int memoryflag,
                   int *sendsize, int *recvsize)
{
  Remap3d *remap = (Remap3d *) ptr;
  remap->setup(in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
               out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
               nqty,permute,memoryflag,*sendsize,*recvsize);
}

/* ----------------------------------------------------------------------
   perform a 3d Remap
------------------------------------------------------------------------- */

void remap3d_remap(void *ptr, FFT_SCALAR *in, FFT_SCALAR *out,
                   FFT_SCALAR *sendbuf, FFT_SCALAR *recvbuf)
{
  Remap3d *remap = (Remap3d *) ptr;
  remap->remap(in,out,sendbuf,recvbuf);
}
