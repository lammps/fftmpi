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

// Compute a forward/inverse complex FFT using fftMPI
//   change FFT size by editing 3 "FFT size" lines
//   run on any number of procs

// Run syntax:
// % simple_c               # run in serial
// % mpirun -np 4 simple_c  # run in parallel

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#include "fft3d_wrap.h"

// FFT size

#define NFAST 128
#define NMID 128
#define NSLOW 128

// precision-dependent settings

#ifdef FFT_SINGLE
int precision = 1;
#else
int precision = 2;
#endif

// main program

int main(int narg, char **args)
{
  // setup MPI

  MPI_Init(&narg,&args);
  MPI_Comm world = MPI_COMM_WORLD;

  int me,nprocs;
  MPI_Comm_size(world,&nprocs);
  MPI_Comm_rank(world,&me);

  // instantiate FFT

  void *fft;
  fft3d_create(world,precision,&fft);

  // simple algorithm to factor Nprocs into roughly cube roots

  int npfast,npmid,npslow;

  npfast = (int) pow(nprocs,1.0/3.0);
  while (npfast < nprocs) {
    if (nprocs % npfast == 0) break;
    npfast++;
  }
  int npmidslow = nprocs / npfast;
  npmid = (int) sqrt(npmidslow);
  while (npmid < npmidslow) {
    if (npmidslow % npmid == 0) break;
    npmid++;
  }
  npslow = nprocs / npfast / npmid;

  // partition grid into Npfast x Npmid x Npslow bricks

  int nfast,nmid,nslow;
  int ilo,ihi,jlo,jhi,klo,khi;

  nfast = NFAST;
  nmid = NMID;
  nslow = NSLOW;

  int ipfast = me % npfast;
  int ipmid = (me/npfast) % npmid;
  int ipslow = me / (npfast*npmid);

  ilo = (int) 1.0*ipfast*nfast/npfast;
  ihi = (int) 1.0*(ipfast+1)*nfast/npfast - 1;
  jlo = (int) 1.0*ipmid*nmid/npmid;
  jhi = (int) 1.0*(ipmid+1)*nmid/npmid - 1;
  klo = (int) 1.0*ipslow*nslow/npslow;
  khi = (int) 1.0*(ipslow+1)*nslow/npslow - 1;

  // setup FFT, could replace with tune()

  int fftsize,sendsize,recvsize;
  fft3d_setup(fft,nfast,nmid,nslow,
              ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi,
              0,&fftsize,&sendsize,&recvsize);

  // tune FFT, could replace with setup()

  //fft3d_tune(fft,nfast,nmid,nslow,
  //           ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi,
  //           0,&fftsize,&sendsize,&recvsize,0,5,10.0,0);

  // initialize each proc's local grid
  // global initialization is specific to proc count

  FFT_SCALAR *work = (FFT_SCALAR *) malloc(2*fftsize*sizeof(FFT_SCALAR));

  int n = 0;
  for (int k = klo; k <= khi; k++) {
    for (int j = jlo; j <= jhi; j++) {
      for (int i = ilo; i <= ihi; i++) {
        work[n] = (double) n;
        n++;
        work[n] = (double) n;
        n++;
      }
    }
  }

  // perform 2 FFTs

  double timestart = MPI_Wtime();
  fft3d_compute(fft,work,work,1);        // forward FFT
  fft3d_compute(fft,work,work,-1);       // inverse FFT
  double timestop = MPI_Wtime();

  if (me == 0) {
    printf("Two %dx%dx%d FFTs on %d procs as %dx%dx%d grid\n",
           nfast,nmid,nslow,nprocs,npfast,npmid,npslow);
    printf("CPU time = %g secs\n",timestop-timestart);
  }

  // find largest difference between initial/final values
  // should be near zero

  n = 0;
  double mydiff = 0.0;
  for (int k = klo; k <= khi; k++) {
    for (int j = jlo; j <= jhi; j++) {
      for (int i = ilo; i <= ihi; i++) {
        if (fabs(work[n]-n) > mydiff) mydiff = fabs(work[n]-n);
        n++;
        if (fabs(work[n]-n) > mydiff) mydiff = fabs(work[n]-n);
        n++;
      }
    }
  }

  double alldiff;
  MPI_Allreduce(&mydiff,&alldiff,1,MPI_DOUBLE,MPI_MAX,world);  
  if (me == 0) printf("Max difference in initial/final values = %g\n",alldiff);

  // clean up

  free(work);
  fft3d_destroy(fft);
  MPI_Finalize();
}
