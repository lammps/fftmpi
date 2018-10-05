! /* ----------------------------------------------------------------------
!    fftMPI - library for computing 3d/2d FFTs in parallel
!    http://fftmpi.sandia.gov, Sandia National Laboratories
!    Steve Plimpton, sjplimp@sandia.gov
!
!    Copyright 2018 National Technology & Engineering Solutions of
!    Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
!    NTESS, the U.S. Government retains certain rights in this software.
!    This software is distributed under the modified Berkeley Software
!    Distribution (BSD) License.
!
!    See the README file in the top-level fftMPI directory.
! ------------------------------------------------------------------------- */

! Compute a forward/inverse complex FFT using fftMPI
!   change FFT size by editing 3 "FFT size" lines
!   run on any number of procs

! Run syntax:
! % simple_f90               # run in serial
! % mpirun -np 4 simple_f90  # run in parallel

! main program

program simple_f90

use iso_c_binding
use fft3d_wrap
implicit none

include 'mpif.h'

! data declarations

integer world,me,nprocs
integer i,j,k,n,precision,ierr
integer nfast_user,nmid_user,nslow_user
integer nfast,nmid,nslow
integer npfast,npmid,npslow,npmidslow,ipfast,ipmid,ipslow
integer ilo,ihi,jlo,jhi,klo,khi
integer fftsize,sendsize,recvsize
real*8 timestart,timestop,mydiff,alldiff
type(c_ptr) :: fft

#ifdef FFT_SINGLE
real(4), allocatable, target :: work(:)
#else
real(8), allocatable, target :: work(:)
#endif

! fft size

nfast_user = 128
nmid_user = 128
nslow_user = 128

! setup mpi

call MPI_Init(ierr)
world = MPI_COMM_WORLD

call MPI_Comm_size(world,nprocs,ierr)
call MPI_Comm_rank(world,me,ierr)

! instantiate fft

#ifdef FFT_SINGLE
precision = 1;
#else
precision = 2;
#endif

call fft3d_create(world,precision,fft)

! simple algorithm to factor nprocs into roughly cube roots

npfast = nprocs**(1.0/3.0)
do while (npfast < nprocs)
  if (mod(nprocs,npfast) == 0) exit
  npfast = npfast + 1
enddo

npmidslow = nprocs / npfast
npmid = sqrt(1.0*npmidslow)
do while (npmid < npmidslow)
  if (mod(npmidslow,npmid) == 0) exit
  npmid = npmid + 1
enddo
npslow = nprocs / npfast / npmid

! partition grid into npfast x npmid x npslow bricks

nfast = nfast_user
nmid = nmid_user
nslow = nslow_user

ipfast = mod(me,npfast)
ipmid = mod((me/npfast),npmid)
ipslow = me / (npfast*npmid)

ilo = 1.0*ipfast*nfast/npfast + 1
ihi = 1.0*(ipfast+1)*nfast/npfast
jlo = 1.0*ipmid*nmid/npmid + 1
jhi = 1.0*(ipmid+1)*nmid/npmid
klo = 1.0*ipslow*nslow/npslow + 1
khi = 1.0*(ipslow+1)*nslow/npslow

! setup fft, could replace with tune()

call fft3d_setup(fft,nfast,nmid,nslow, &
        ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi, &
        0,fftsize,sendsize,recvsize)

! tune fft, could replace with setup()

!call fft3d_tune(fft,nfast,nmid,nslow, &
!        ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi, &
!        0,fftsize,sendsize,recvsize,0,5,10.0d0,0)

! initialize each proc's local grid
! global initialization is specific to proc count

allocate(work(2*fftsize))

n = 1
do k = klo,khi
  do j = jlo,jhi
    do i = ilo,ihi
      work(n) = n
      n = n + 1
      work(n) = n
      n = n + 1
    enddo
  enddo
enddo

! perform 2 ffts

timestart = mpi_wtime()
call fft3d_compute(fft,c_loc(work),c_loc(work),1)        ! forward fft
call fft3d_compute(fft,c_loc(work),c_loc(work),-1)       ! inverse fft
timestop = mpi_wtime()

if (me == 0) then
  print *,"Two",nfast,"x",nmid,"x",nslow,"ffts on",nprocs, &
          "procs as ",npfast,"x",npmid,"x",npslow,"grid"
  print *,"CPU time =",timestop-timestart,"secs"
endif

! find largest difference between initial/final values
! should be near zero

n = 1
mydiff = 0.0
do k = klo,khi
  do j = jlo,jhi
    do i = ilo,ihi
      if (abs(work(n)-n) > mydiff) mydiff = abs(work(n)-n)
      n = n + 1
      if (abs(work(n)-n) > mydiff) mydiff = abs(work(n)-n)
      n = n + 1
    enddo
  enddo
enddo

call MPI_Allreduce(mydiff,alldiff,1,mpi_double,mpi_max,world,ierr)  
if (me == 0) print *,"Max difference in initial/final values =",alldiff

! clean up

deallocate(work)
call fft3d_destroy(fft)
call mpi_finalize()

end program simple_f90
