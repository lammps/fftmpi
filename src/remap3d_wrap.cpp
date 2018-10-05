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

// C interface to fftMPI library, 3d Remap functions

#include <string.h>
#include <stdlib.h>

#include "remap3d_wrap.h"
#include "remap3d.h"

using namespace FFTMPI_NS;

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
   create plan for performing a 3d Remap
   Fortran interface where indices are 1 to N inclusive
------------------------------------------------------------------------- */

void remap3d_setup_fortran(void *ptr,
                           int in_ilo, int in_ihi, int in_jlo, 
                           int in_jhi, int in_klo, int in_khi,
                           int out_ilo, int out_ihi, int out_jlo, 
                           int out_jhi, int out_klo, int out_khi,
                           int nqty, int permute, int memoryflag,
                           int *sendsize, int *recvsize)
{
  Remap3d *remap = (Remap3d *) ptr;
  remap->setup(in_ilo-1,in_ihi-1,in_jlo-1,in_jhi-1,in_klo-1,in_khi-1,
               out_ilo-1,out_ihi-1,out_jlo-1,out_jhi-1,out_klo-1,out_khi-1,
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
