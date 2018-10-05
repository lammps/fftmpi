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

// test driver on 3d FFT from FFT library
// for benchmarking, timing purposes

// Syntax: test3d switch args switch args ...
//      -h = print help message
//      -g Nx Ny Nz = grid size (default = 8 8 8)
//      -pin Px Py Pz = proc grid (default = 0 0 0)
//         specify 3d grid of procs for initial partition
//         0 0 0 = code chooses Px Py Pz, will be bricks
//      -pout Px Py Pz = proc grid (default = 0 0 0)
//         specify 3d grid of procs for final partition
//         0 0 0 = code chooses Px Py Pz
//                 will be bricks for mode = 0/2
//                 will be z pencils for mode = 1/3
//      -n Nloop = iteration count (default = 1)
//                 can use 0 if -tune enabled, then will set by tuning operation
//      -m 0/1/2/3 = FFT mode (default = 0)
//         0 = 1 iteration = forward full FFT, backward full FFT
//         1 = 1 iteration = forward convolution FFT, backward convolution FFT
//         2 = 1 iteration = just forward full FFT
//         3 = 1 iteration = just forward convolution FFT
//         full FFT returns data to original layout
//         forward convolution FFT is brick -> z-pencil
//         backward convolution FFT is z-pencil -> brick
//      -i zero/step/index/82783 = initialization (default = zero)
//         zero = initialize grid with 0.0
//         step = initialize with 3d step function
//         index = ascending integers 1 to Nx+Ny+Nz
//         digits = random number seed
//      -tune nper tmax extra
//         nper = # of FFTs per trial run
//         tmax = tune within tmax CPU secs, 0.0 = unlimited
//         extra = 1 for detailed timing of trial runs, else 0
//      -c point/all/combo = communication flag (default = point)
//         point = point-to-point comm
//         all = use MPI_all2all collective
//         combo = point for pencil2brick, all2all for pencil2pencil
//      -e pencil/brick = exchange flag (default = pencil)
//         pencil = pencil to pencil data exchange (4 stages for full FFT)
//         brick = brick to pencil data exchange (6 stages for full FFT)
//      -p array/ptr/memcpy
//         pack/unpack methods for data remapping (default = memcpy)
//         array = array based
//         ptr = pointer based
//         memcpy = memcpy based
//      -t = provide more timing details (not set by default)
//         include timing breakdown, not just summary
//      -r = remap only, no 1d FFTs (not set by default)
//         useful for debugging
//      -o = output initial/final grid (not set by default)
//         only useful for small problems
//      -v = verify correctness of answer (not set by default)
//         only possible for FFT mode = 0/1

// include files

#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#include "fft3d.h"

using namespace FFTMPI_NS;

// precision-dependent settings

#ifdef FFT_SINGLE
int precision = 1;
#else
int precision = 2;
#endif

// memory alignment settings

#if defined(__INTEL_COMPILER)
#ifndef FFT_INTEL_NO_TBB
#define FFT_USE_TBB_ALLOCATOR
#include "tbb/scalable_allocator.h"
#else
#include <malloc.h>
#endif
#endif

#if !defined(FFT_MEMALIGN)
#define FFT_MEMALIGN 64
#endif

// RNG constants

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

// common data

typedef int64_t bigint;

int me,nprocs;
MPI_Comm world;

int nx,ny,nz;
int inpx,inpy,inpz,outpx,outpy,outpz;
int nloop;
int mode;
int iflag,oflag,tflag,cflag,eflag,pflag,rflag,vflag;
int seed,seedinit;
int tuneflag,tuneper,tuneextra;
double tunemax;

int inxlo,inxhi,inylo,inyhi,inzlo,inzhi;         // initial partition of grid
int outxlo,outxhi,outylo,outyhi,outzlo,outzhi;   // final partition of grid
int nfft_in;                      // # of grid pts I own in initial partition
int nfft_out;                     // # of grid pts I own in final partition
int fftsize;                      // FFT buffer size returned by FFT setup

FFT3d *fft;
double timefft,timeinit,timesetup,timetune;
double epsmax;

FFT_SCALAR *work;

// functions

void options(int, char **);
void proc_setup(int);
void proc3d(int &, int &, int &);
void proc2d(int &, int &, int &);
void grid_setup();
void plan();
void allocate();
void initialize();
void output(int, const char *);
void validate();
void timing();
void deallocate();

double random(int &);

void error_all(const char *);
void error_one(const char *);

void *smalloc(int64_t);
void sfree(void *);

// constants

const char *syntax = 
  "Syntax: test3d -g Nx Nx Nz -p Px Py Pz -n Nloop -m 0/1/2/3\n"
  "               -i zero/step/82783 -m 0/1/2/3 -tune nper tmax extra\n"
  "               -c point/all/combo -e pencil/brick -p array/ptr/memcpy\n"
  "               -t -r -o -v";

enum{ZERO,STEP,INDEX,RANDOM};
enum{POINT,ALL2ALL,COMBO};
enum{PENCIL,BRICK};
enum{ARRAY,POINTER,MEMCPY};
enum{IN,OUT};

/* ----------------------------------------------------------------------
   main progam
------------------------------------------------------------------------- */

int main(int narg, char **args)
{
  // MPI setup

  MPI_Init(&narg,&args);
  world = MPI_COMM_WORLD;

  MPI_Comm_size(world,&nprocs);
  MPI_Comm_rank(world,&me);

  // parse command-line args

  options(narg,args);

  // partition FFT grid across procs, for both input and output
  // create FFT plan, tune if requested
  // allocate grid
  // initialize FFT grid
  // grid output

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  proc_setup(IN);
  proc_setup(OUT);
  grid_setup();
  plan();
  allocate();
  initialize();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();
  timeinit = time2 - time1;

  if (oflag) output(0,"Initial grid");

  // perform FFTs

  MPI_Barrier(world);
  time1 = MPI_Wtime();

  if (mode < 2) {
    for (int i = 0; i < nloop; i++) {
      fft->compute(work,work,1);
      //if (oflag) output(1,"Middle grid");
      fft->compute(work,work,-1);
    }
  } else {
    for (int i = 0; i < nloop; i++) {
      fft->compute(work,work,1);
    }
  }
  
  MPI_Barrier(world);
  time2 = MPI_Wtime();
  timefft = time2 - time1;

  // validation check on result
  // grid output
  // timing results
  // deallocate grid and plan

  if (vflag) validate();
  if (oflag) {
    if (mode < 2) output(0,"Final grid");
    else output(1,"Final grid");
  }
  timing();
  deallocate();
  delete fft;

  // shut down MPI

  MPI_Finalize();
}

/* ----------------------------------------------------------------------
   parse command-line options
   all options have defaults
------------------------------------------------------------------------- */

void options(int narg, char **args)
{
  // defaults

  nx = ny = nz = 8;
  inpx = inpy = inpz = 0;
  outpx = outpy = outpz = 0;
  nloop = 1;
  iflag = ZERO;
  tuneflag = 0;
  mode = 0;
  cflag = COMBO;
  eflag = PENCIL;
  pflag = MEMCPY;
  tflag = 0;
  rflag = 0;
  oflag = 0;
  vflag = 0;

  // parse args

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(args[iarg],"-h") == 0) {
      error_all(syntax);
    } else if (strcmp(args[iarg],"-g") == 0) {
      if (iarg+4 > narg) error_all(syntax);
      nx = atoi(args[iarg+1]);
      ny = atoi(args[iarg+2]);
      nz = atoi(args[iarg+3]);
      iarg += 4;
    } else if (strcmp(args[iarg],"-pin") == 0) {
      if (iarg+4 > narg) error_all(syntax);
      inpx = atoi(args[iarg+1]);
      inpy = atoi(args[iarg+2]);
      inpz = atoi(args[iarg+3]);
      iarg += 4;
    } else if (strcmp(args[iarg],"-pout") == 0) {
      if (iarg+4 > narg) error_all(syntax);
      outpx = atoi(args[iarg+1]);
      outpy = atoi(args[iarg+2]);
      outpz = atoi(args[iarg+3]);
      iarg += 4;
    } else if (strcmp(args[iarg],"-n") == 0) {
      if (iarg+2 > narg) error_all(syntax);
      nloop = atoi(args[iarg+1]);
      iarg += 2;
    } else if (strcmp(args[iarg],"-i") == 0) {
      if (iarg+2 > narg) error_all(syntax);
      if (strcmp(args[iarg+1],"zero") == 0) iflag = ZERO;
      else if (strcmp(args[iarg+1],"step") == 0) iflag = STEP;
      else if (strcmp(args[iarg+1],"index") == 0) iflag = INDEX;
      else {
        iflag = RANDOM;
        // per-processor RNG seed
        seed = seedinit = atoi(args[iarg+1]) + me;
      }
      iarg += 2;
    } else if (strcmp(args[iarg],"-tune") == 0) {
      if (iarg+4 > narg) error_all(syntax);
      tuneflag = 1;
      tuneper = atoi(args[iarg+1]);
      tunemax = atof(args[iarg+2]);
      tuneextra = atoi(args[iarg+3]);
      iarg += 4;
    } else if (strcmp(args[iarg],"-m") == 0) {
      if (iarg+2 > narg) error_all(syntax);
      mode = atoi(args[iarg+1]);
      iarg += 2;
    } else if (strcmp(args[iarg],"-c") == 0) {
      if (iarg+2 > narg) error_all(syntax);
      if (strcmp(args[iarg+1],"point") == 0) cflag = POINT;
      else if (strcmp(args[iarg+1],"all") == 0) cflag = ALL2ALL;
      else if (strcmp(args[iarg+1],"combo") == 0) cflag = COMBO;
      else error_all(syntax);
      iarg += 2;
    } else if (strcmp(args[iarg],"-e") == 0) {
      if (iarg+2 > narg) error_all(syntax);
      if (strcmp(args[iarg+1],"pencil") == 0) eflag = PENCIL;
      else if (strcmp(args[iarg+1],"brick") == 0) eflag = BRICK;
      else error_all(syntax);
      iarg += 2;
    } else if (strcmp(args[iarg],"-p") == 0) {
      if (iarg+2 > narg) error_all(syntax);
      if (strcmp(args[iarg+1],"array") == 0) pflag = ARRAY;
      else if (strcmp(args[iarg+1],"ptr") == 0) pflag = POINTER;
      else if (strcmp(args[iarg+1],"memcpy") == 0) pflag = MEMCPY;
      else error_all(syntax);
      iarg += 2;
    } else if (strcmp(args[iarg],"-t") == 0) {
      tflag = 1;
      iarg += 1;
    } else if (strcmp(args[iarg],"-r") == 0) {
      rflag = 1;
      iarg += 1;
    } else if (strcmp(args[iarg],"-o") == 0) {
      oflag = 1;
      iarg += 1;
    } else if (strcmp(args[iarg],"-v") == 0) {
      vflag = 1;
      iarg += 1;
    } else error_all(syntax);
  }

  // sanity check on args

  if (nx <= 0 || ny <= 0 || nz <= 0) error_all("Invalid grid size");

  if (inpx == 0 && inpy == 0 && inpz == 0);
  else if (inpx <= 0 || inpy <= 0 || inpz <= 0) error_all("Invalid proc grid");
  else if (inpx*inpy*inpz != nprocs) 
    error_all("Specified proc grid does not match nprocs");

  if (outpx == 0 && outpy == 0 && outpz == 0);
  else if (outpx <= 0 || outpy <= 0 || outpz <= 0) 
    error_all("Invalid proc grid");
  else if (outpx*outpy*outpz != nprocs) 
    error_all("Specified proc grid does not match nprocs");

  if (nloop < 0) error_all("Invalid Nloop");
  if (nloop == 0 && tuneflag == 0) error_all("Invalid Nloop");
  if (iflag == RANDOM && seed <= 0) error_all("Invalid initialize setting");
  if (mode < 0 || mode > 3) error_all("Invalid FFT mode");
  if (mode > 1 && vflag) error_all("Cannot validate forward only FFT");

  if (tuneflag && tuneper <= 0) error_all("Invalid tune nper");
  if (tuneflag && tunemax < 0.0) error_all("Invalid tune tmax");
  if (tuneflag && (tuneextra < 0 || tuneextra > 1))
    error_all("Invalid tune extra");
  if (tuneflag && rflag) error_all("Cannot tune with remap only");
}

/* ----------------------------------------------------------------------
   partition processors across grid dimensions
   flag = IN for input partitions, or OUT for output partitions
   if user set Px,Py,Pz -> just return
   for IN:
     assign nprocs as bricks to 3d grid to minimize surface area per proc
     derived from SPPARKS Domain::procs2domain_3d()
   for OUT:
     assign nprocs as rectangles to xy grid to minimize surface area per proc
     derived from SPPARKS Domain::procs2domain_2d()
------------------------------------------------------------------------- */

void proc_setup(int flag)
{
  if (flag == IN) {
    if (inpx != 0 || inpy != 0 || inpz != 0) return;
    proc3d(inpx,inpy,inpz);
  }

  if (flag == OUT) {
    if (outpx != 0 || outpy != 0 || outpz != 0) return;
    if (mode == 0 || mode == 2) proc3d(outpx,outpy,outpz);
    if (mode == 1 || mode == 3) proc2d(outpx,outpy,outpz);
  }
}

void proc3d(int &px, int &py, int &pz)
{
  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;
  double xprd = nx;
  double yprd = ny;
  double zprd = nz;
  
  double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);
  
  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  
  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      nremain = nprocs/ipx;
      ipy = 1;
      while (ipy <= nremain) {
        if (nremain % ipy == 0) {
          ipz = nremain/ipy;
          boxx = xprd/ipx;
          boxy = yprd/ipy;
          boxz = zprd/ipz;
          surf = boxx*boxy + boxy*boxz + boxz*boxx;
          if (surf < bestsurf) {
            bestsurf = surf;
            px = ipx;
            py = ipy;
            pz = ipz;
          }
        }
        ipy++;
      }
    }
    ipx++;
  }

  if (px*py*pz != nprocs) 
    error_all("Computed proc grid does not match nprocs");
}

void proc2d(int &px, int &py, int &pz)
{
  int ipx,ipy;
  double boxx,boxy,surf;
  double xprd = nx;
  double yprd = ny;
  
  double bestsurf = 2.0 * (xprd+yprd);
  
  // loop thru all possible factorizations of nprocs
  // surf = surface area of a proc sub-domain
  
  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = xprd/ipx;
      boxy = yprd/ipy;
      surf = boxx + boxy;
      if (surf < bestsurf) {
        bestsurf = surf;
        px = ipx;
        py = ipy;
      }
    }
    ipx++;
  }
  
  pz = 1;
  if (px*py*pz != nprocs) 
    error_all("Computed proc grid does not match nprocs");
}

/* ----------------------------------------------------------------------
   partition FFT grid
   once for input grid, once for output grid
   use Px,Py,Pz for in/out
------------------------------------------------------------------------- */

void grid_setup()
{
  // ipx,ipy,ipz = my position in input 3d grid of procs

  int ipx = me % inpx;
  int ipy = (me/inpx) % inpy;
  int ipz = me / (inpx*inpy);

  // nlo,nhi = lower/upper limits of the 3d brick I own

  inxlo = static_cast<int> (1.0 * ipx * nx / inpx);
  inxhi = static_cast<int> (1.0 * (ipx+1) * nx / inpx) - 1;

  inylo = static_cast<int> (1.0 * ipy * ny / inpy);
  inyhi = static_cast<int> (1.0 * (ipy+1) * ny / inpy) - 1;

  inzlo = static_cast<int> (1.0 * ipz * nz / inpz);
  inzhi = static_cast<int> (1.0 * (ipz+1) * nz / inpz) - 1;

  nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1);

  // ipx,ipy,ipz = my position in output 3d grid of procs

  ipx = me % outpx;
  ipy = (me/outpx) % outpy;
  ipz = me / (outpx*outpy);

  // nlo,nhi = lower/upper limits of the 3d brick I own

  outxlo = static_cast<int> (1.0 * ipx * nx / outpx);
  outxhi = static_cast<int> (1.0 * (ipx+1) * nx / outpx) - 1;

  outylo = static_cast<int> (1.0 * ipy * ny / outpy);
  outyhi = static_cast<int> (1.0 * (ipy+1) * ny / outpy) - 1;

  outzlo = static_cast<int> (1.0 * ipz * nz / outpz);
  outzhi = static_cast<int> (1.0 * (ipz+1) * nz / outpz) - 1;

  nfft_out = (outxhi-outxlo+1) * (outyhi-outylo+1) * (outzhi-outzlo+1);
}

/* ----------------------------------------------------------------------
   create FFT plan
------------------------------------------------------------------------- */

void plan()
{
  fft = new FFT3d(world,precision);
  fft->remaponly = rflag;

  fft->collective = cflag;
  fft->exchange = eflag;
  fft->packflag = pflag;

  int permute;
  if (mode == 0 || mode == 2) permute = 0;
  else permute = 2;

  // will use fftsize to allocate work buffer
  // ignore sendsize, recvsize b/c let FFT allocate remap buffers internally
  // set timesetup and timetune
  // reset nloop if tuning and user nloop = 0

  int sendsize,recvsize;

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  if (!tuneflag) {
    fft->setup(nx,ny,nz,
               inxlo,inxhi,inylo,inyhi,inzlo,inzhi,
               outxlo,outxhi,outylo,outyhi,outzlo,outzhi,
               permute,fftsize,sendsize,recvsize);
  } else {
    int flag = 0;
    if (mode >= 2) flag = 1;
    fft->tune(nx,ny,nz,
              inxlo,inxhi,inylo,inyhi,inzlo,inzhi,
              outxlo,outxhi,outylo,outyhi,outzlo,outzhi,
              permute,fftsize,sendsize,recvsize,
              flag,tuneper,tunemax,tuneextra);
    if (nloop == 0) nloop = fft->npertrial;
  }

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  if (!tuneflag) {
    timesetup = time2 - time1;
    timetune = 0.0;
  } else {
    timesetup = fft->setuptime;
    timetune = time2 - time1;
  }
}

/* ----------------------------------------------------------------------
   allocate memory for FFT grid
------------------------------------------------------------------------- */

void allocate()
{
  bigint nbytes = ((bigint) sizeof(FFT_SCALAR)) * 2*fftsize;
  work = (FFT_SCALAR *) smalloc(nbytes);
  if (nbytes && work == NULL) error_one("Failed malloc for FFT grid");
}

/* ----------------------------------------------------------------------
   initialize FFT grid
------------------------------------------------------------------------- */

void initialize()
{
  if (iflag == ZERO) {
    for (int m = 0; m < 2*nfft_in; m++) work[m] = 0.0;

  } else if (iflag == STEP) {
    int ilocal,jlocal,klocal,iglobal,jglobal,kglobal;
    int nxlocal = inxhi - inxlo + 1;
    int nylocal = inyhi - inylo + 1;

    for (int m = 0; m < nfft_in; m++) {
      ilocal = m % nxlocal;
      jlocal = (m/nxlocal) % nylocal;
      klocal = m / (nxlocal*nylocal);
      iglobal = inxlo + ilocal;
      jglobal = inylo + jlocal;
      kglobal = inzlo + klocal;
      if (iglobal < nx/2 && jglobal < ny/2 && kglobal < nz/2) work[2*m] = 1.0;
      else work[2*m] = 0.0;
      work[2*m+1] = 0.0;
    }

  } else if (iflag == INDEX) {
    int ilocal,jlocal,klocal,iglobal,jglobal,kglobal;
    int nxlocal = inxhi - inxlo + 1;
    int nylocal = inyhi - inylo + 1;

    for (int m = 0; m < nfft_in; m++) {
      ilocal = m % nxlocal;
      jlocal = (m/nxlocal) % nylocal;
      klocal = m / (nxlocal*nylocal);
      iglobal = inxlo + ilocal;
      jglobal = inylo + jlocal;
      kglobal = inzlo + klocal;
      work[2*m] = kglobal + jglobal + iglobal + 1;
      work[2*m+1] = 0.0;
    }

  } else if (iflag == RANDOM) {
    for (int m = 0; m < 2*nfft_in; m++) 
      work[m] = random(seed);
  }
}

/* ----------------------------------------------------------------------
   output FFT grid values
   flag = 0 for initial partition
   flag = 1 for final partition
------------------------------------------------------------------------- */

void output(int flag, const char *str)
{
  int tmp;

  if (me == 0) printf("%s\n",str);

  for (int iproc = 0; iproc < nprocs; iproc++) {
    if (me != iproc) continue;
    if (me >= 1) MPI_Recv(&tmp,0,MPI_INT,me-1,0,world,MPI_STATUS_IGNORE);

    int ilocal,jlocal,klocal,iglobal,jglobal,kglobal;

    if (flag == 0) {
      int nxlocal = inxhi - inxlo + 1;
      int nylocal = inyhi - inylo + 1;

      for (int m = 0; m < nfft_in; m++) {
        ilocal = m % nxlocal;
        jlocal = (m/nxlocal) % nylocal;
        klocal = m / (nxlocal*nylocal);
        iglobal = inxlo + ilocal;
        jglobal = inylo + jlocal;
        kglobal = inzlo + klocal;
        printf("Value (%d,%d,%d) on proc %d = (%g,%g)\n",
               iglobal,jglobal,kglobal,
               me,work[2*m],work[2*m+1]);
      } 
    } else {
      int nxlocal = outxhi - outxlo + 1;
      int nylocal = outyhi - outylo + 1;

      for (int m = 0; m < nfft_in; m++) {
        ilocal = m % nxlocal;
        jlocal = (m/nxlocal) % nylocal;
        klocal = m / (nxlocal*nylocal);
        iglobal = outxlo + ilocal;
        jglobal = outylo + jlocal;
        kglobal = outzlo + klocal;
        printf("Value (%d,%d,%d) on proc %d = (%g,%g)\n",
               iglobal,jglobal,kglobal,
               me,work[2*m],work[2*m+1]);
      }
    }

    if (me < nprocs-1) MPI_Send(&tmp,0,MPI_INT,me+1,0,world);
  }
}

/* ----------------------------------------------------------------------
   validation check for correct result
------------------------------------------------------------------------- */

void validate()
{
  double delta;
  double epsilon = 0.0;

  if (iflag == ZERO) {
    for (int m = 0; m < 2*nfft_in; m++) {
      delta = fabs(work[m]);
      if (delta > epsilon) epsilon = delta;
    }

  } else if (iflag == STEP) {
    int ilocal,jlocal,klocal,iglobal,jglobal,kglobal;
    int nxlocal = inxhi - inxlo + 1;
    int nylocal = inyhi - inylo + 1;
    double value;

    for (int m = 0; m < nfft_in; m++) {
      ilocal = m % nxlocal;
      jlocal = (m/nxlocal) % nylocal;
      klocal = m / (nxlocal*nylocal);
      iglobal = inxlo + ilocal;
      jglobal = inylo + jlocal;
      kglobal = inzlo + klocal;
      if (iglobal < nx/2 && jglobal < ny/2 && kglobal < nz/2) value = 1.0;
      else value = 0.0;
      delta = fabs(work[2*m]-value);
      if (delta > epsilon) epsilon = delta;
      delta = fabs(work[2*m+1]);
      if (delta > epsilon) epsilon = delta;
    }

  } else if (iflag == INDEX) {
    int ilocal,jlocal,klocal,iglobal,jglobal,kglobal;
    int nxlocal = inxhi - inxlo + 1;
    int nylocal = inyhi - inylo + 1;
    double value;

    for (int m = 0; m < nfft_in; m++) {
      ilocal = m % nxlocal;
      jlocal = (m/nxlocal) % nylocal;
      klocal = m / (nxlocal*nylocal);
      iglobal = inxlo + ilocal;
      jglobal = inylo + jlocal;
      kglobal = inzlo + klocal;
      value = kglobal + jglobal + iglobal + 1;
      delta = fabs(work[2*m]-value);
      if (delta > epsilon) epsilon = delta;
      delta = fabs(work[2*m+1]);
      if (delta > epsilon) epsilon = delta;
    }

  } else if (iflag == RANDOM) {
    double newvalue;
    seed = seedinit;
    for (int m = 0; m < 2*nfft_in; m++) {
      newvalue = random(seed);
      delta = fabs(work[m]-newvalue);
      if (delta > epsilon) epsilon = delta;
    }
  }

  MPI_Allreduce(&epsilon,&epsmax,1,MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   output timing data
------------------------------------------------------------------------- */

void timing()
{
  double time1d,time_remap;
  double time_remap1,time_remap2,time_remap3,time_remap4;

  // perform only 1d FFTs

  if (tflag) {
    for (int i = 0; i < 2*nfft_in; i++) work[i] = 0.0;

    MPI_Barrier(world);
    double time1 = MPI_Wtime();

    if (mode < 2) {
      for (int i = 0; i < nloop; i++) {
        fft->only_1d_ffts(work,1);
        fft->only_1d_ffts(work,-1);
      }
    } else {
      for (int i = 0; i < nloop; i++) {
        fft->only_1d_ffts(work,1);
      }
    }

    MPI_Barrier(world);
    double time2 = MPI_Wtime();
    time1d = time2 - time1;
  }

  // perform all remaps

  if (tflag) {
    for (int i = 0; i < 2*nfft_in; i++) work[i] = 0.0;

    MPI_Barrier(world);
    double time1 = MPI_Wtime();

    if (mode < 2) {
      for (int i = 0; i < nloop; i++) {
        fft->only_remaps(work,work,1);
        fft->only_remaps(work,work,-1);
      }
    } else {
      for (int i = 0; i < nloop; i++) {
        fft->only_remaps(work,work,1);
      }
    }

    MPI_Barrier(world);
    double time2 = MPI_Wtime();
    time_remap = time2 - time1;
  }

  // perform only single remaps

  if (tflag) {
    for (int i = 0; i < 2*nfft_in; i++) work[i] = 0.0;

    MPI_Barrier(world);
    double time1 = MPI_Wtime();

    if (mode < 2) {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,1);
        fft->only_one_remap(work,work,-1,1);
      }
    } else {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,1);
      }
    }

    MPI_Barrier(world);
    double time2 = MPI_Wtime();
    time_remap1 = time2 - time1;

    if (mode < 2) {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,2);
        fft->only_one_remap(work,work,-1,2);
      }
    } else {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,2);
      }
    }

    MPI_Barrier(world);
    double time3 = MPI_Wtime();
    time_remap2 = time3 - time2;

    if (mode < 2) {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,3);
        fft->only_one_remap(work,work,-1,3);
      }
    } else {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,3);
      }
    }

    MPI_Barrier(world);
    double time4 = MPI_Wtime();
    time_remap3 = time4 - time3;

    if (mode < 2) {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,4);
        fft->only_one_remap(work,work,-1,4);
      }
    } else {
      for (int i = 0; i < nloop; i++) {
        fft->only_one_remap(work,work,1,4);
      }
    }

    MPI_Barrier(world);
    double time5 = MPI_Wtime();
    time_remap4 = time5 - time4;
  }

  // stats output
  // nfft = # of FFTs performed = 2x larger for modes 0,1

  int nfft;
  if (mode < 2) nfft = 2*nloop;
  else nfft = nloop;

  double onetime = timefft/nfft;
  double nsize = 1.0 * nx * ny * nz;
  double log2n = log(nsize)/log(2.0);
  double floprate = 5.0 * nsize * log2n / onetime / (1024*1024*1024);
  bigint gridbytes = ((bigint) sizeof(FFT_SCALAR)) * 2*fftsize;

  if (me == 0) {
    printf("3d FFTs with %s library, precision = %s\n",
           fft->fft1d,fft->precision);
    printf("Grid size: %d %d %d\n",nx,ny,nz);
    printf("  initial proc grid: %d %d %d\n",inpx,inpy,inpz);
    printf("  x pencil proc grid: %d %d %d\n",
           fft->npfast1,fft->npfast2,fft->npfast3);
    printf("  y pencil proc grid: %d %d %d\n",
           fft->npmid1,fft->npmid2,fft->npmid3);
    printf("  z pencil proc grid: %d %d %d\n",
           fft->npslow1,fft->npslow2,fft->npslow3);
    printf("  3d brick proc grid: %d %d %d\n",
           fft->npbrick1,fft->npbrick2,fft->npbrick3);
    printf("  final proc grid: %d %d %d\n",outpx,outpy,outpz);

    if (tuneflag) {
      printf("Tuning trials & iterations: %d %d\n",fft->ntrial,fft->npertrial);
      for (int i = 0; i < fft->ntrial; i++)
        printf("  coll exch pack 3dFFT 1dFFT remap r1 r2 r3 r4: "
               "%d %d %d %g %g %g %g %g %g %g\n",
               fft->cflags[i],fft->eflags[i],fft->pflags[i],
               fft->tfft[i],fft->t1d[i],fft->tremap[i],
               fft->tremap1[i],fft->tremap2[i],
               fft->tremap3[i],fft->tremap4[i]);
    }

    if (mode == 0)
      printf("%d forward and %d back FFTs on %d procs\n",nloop,nloop,nprocs);
    else if (mode == 1)
      printf("%d forward and %d back convolution FFTs on %d procs\n",
             nloop,nloop,nprocs);
    else if (mode == 2)
      printf("%d forward FFTs on %d procs\n",nloop,nprocs);
    else if (mode == 3)
      printf("%d forward convolution FFTs on %d procs\n",nloop,nprocs);

    printf("Collective, exchange, pack methods: %d %d %d\n",
           fft->collective,fft->exchange,fft->packflag);
    printf("Memory usage (per-proc) for FFT grid = %g MBytes\n",
           (double) gridbytes / 1024/1024);
    printf("Memory usage (per-proc) by fftMPI = %g MBytes\n",
           (double) fft->memusage / 1024/1024);

    if (vflag) printf("Max error = %g\n",epsmax);
    if (!tuneflag) printf("Initialize grid = %g secs\n",timeinit-timesetup);
    else printf("Initialize grid = %g secs\n",timeinit-timetune);
    printf("FFT setup = %g secs\n",timesetup);
    printf("FFT tune = %g secs\n",timetune);
    printf("Time for 3d FFTs = %g secs\n",timefft);
    printf("  time/fft3d = %g secs\n",onetime);
    printf("  flop rate for 3d FFTs = %g Gflops\n",floprate);
    if (tflag) {
      printf("Time for 1d FFTs only = %g secs\n",time1d);
      printf("  time/fft1d = %g secs\n",time1d/nfft);
      printf("  fraction of time in 1d FFTs = %g\n",time1d/timefft);
    }
    if (tflag) {
      printf("Time for remaps only = %g secs\n",time_remap);
      printf("  fraction of time in remaps = %g\n",time_remap/timefft);
      printf("Time for remap #1 = %g secs\n",time_remap1);
      printf("  fraction of time in remap #1 = %g\n",time_remap1/timefft);
      printf("Time for remap #2 = %g secs\n",time_remap2);
      printf("  fraction of time in remap #2 = %g\n",time_remap2/timefft);
      printf("Time for remap #3 = %g secs\n",time_remap3);
      printf("  fraction of time in remap #3 = %g\n",time_remap3/timefft);
      printf("Time for remap #4 = %g secs\n",time_remap4);
      printf("  fraction of time in remap #4 = %g\n",time_remap4/timefft);
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate memory for FFT grid
------------------------------------------------------------------------- */

void deallocate()
{
  sfree(work);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// utility functions
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

/* ----------------------------------------------------------------------
   simple Park RNG
   pass in non-zero seed
------------------------------------------------------------------------- */

double random(int &seed)
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   must be called by all procs in world
   shuts down MPI and exits
------------------------------------------------------------------------- */

void error_all(const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me == 0) printf("ERROR: %s\n",str);
  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   forces abort of entire world if any proc in world calls
------------------------------------------------------------------------- */

void error_one(const char *str)
{
  int me;
  MPI_Comm_rank(world,&me);
  printf("ERROR on proc %d: %s\n",me,str);
  MPI_Abort(world,1);
}

/* ----------------------------------------------------------------------
   safe malloc
------------------------------------------------------------------------- */

void *smalloc(int64_t nbytes)
{
  if (nbytes == 0) return NULL;

#if defined(FFT_MEMALIGN)
  void *ptr;

#if defined(FFT_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_malloc(nbytes,FFT_MEMALIGN);
#else
  int retval = posix_memalign(&ptr,FFT_MEMALIGN,nbytes);
  if (retval) ptr = NULL;
#endif

#else
  void *ptr = malloc(nbytes);
#endif

  return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void sfree(void *ptr)
{
  if (ptr == NULL) return;

#if defined(FFT_USE_TBB_ALLOCATOR)
  scalable_aligned_free(ptr);
#else
  free(ptr);
#endif

  ptr = NULL;
}
