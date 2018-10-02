! Compute a forward/inverse double precision complex FFT using fftMPI
!   change FFT size by editing 3 "FFT size" lines
!   run on any number of procs

! Run syntax:
! % simple_f90               # run in serial
! % mpirun -np 4 simple_f90  # run in parallel

! FFT size

program simple_f90

use iso_c_binding
use fft3d_wrap
implicit none

include 'mpif.h'

! data declarations

INTEGER world,me,nprocs
INTEGER i,j,k,n,PRECISION,ierr
INTEGER NFAST_USER,NMID_USER,NSLOW_USER
integer nfast,nmid,nslow
INTEGER npfast,npmid,npslow,npmidslow,ipfast,ipmid,ipslow
integer ilo,ihi,jlo,jhi,klo,khi
integer fftsize,sendsize,recvsize
REAL*8 timestart,timestop,mydiff,alldiff
REAL(8), ALLOCATABLE, target :: work(:)
TYPE(C_ptr) :: fft

! FFT size

NFAST_USER = 128
NMID_USER = 128
NSLOW_USER = 128

! setup MPI

call MPI_Init(ierr)
world = MPI_COMM_WORLD

call MPI_Comm_size(world,nprocs,ierr)
call MPI_Comm_rank(world,me,ierr)

! instantiate FFT

PRECISION = 2
call fft3d_create(world,PRECISION,fft)

! simple algorithm to factor Nprocs into roughly cube roots

npfast = nprocs**(1.0/3.0)
do while (npfast < nprocs)
  IF (MOD(nprocs,npfast) == 0) exit
  npfast = npfast + 1
enddo

npmidslow = nprocs / npfast
npmid = SQRT(1.0*npmidslow)
do WHILE (npmid < npmidslow)
  IF (mod(npmidslow,npmid) == 0) exit
  npmid = npmid + 1
enddo
npslow = nprocs / npfast / npmid

! partition grid into Npfast x Npmid x Npslow bricks

nfast = NFAST_USER
nmid = NMID_USER
nslow = NSLOW_USER

ipfast = mod(me,npfast)
ipmid = MOD((me/npfast),npmid)
ipslow = me / (npfast*npmid)

ilo = 1.0*ipfast*nfast/npfast + 1
ihi = 1.0*(ipfast+1)*nfast/npfast
jlo = 1.0*ipmid*nmid/npmid + 1
jhi = 1.0*(ipmid+1)*nmid/npmid
klo = 1.0*ipslow*nslow/npslow + 1
khi = 1.0*(ipslow+1)*nslow/npslow

! setup FFT, could replace with tune()

call fft3d_setup(fft,nfast,nmid,nslow, &
        ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi, &
        0,fftsize,sendsize,recvsize)

! tune FFT, could replace with setup()

!call fft3d_tune(fft,nfast,nmid,nslow, &
!        ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi, &
!        0,fftsize,sendsize,recvsize,0,5,10.0d0,0)

! initialize each proc's local grid
! global initialization is specific to proc count

ALLOCATE(work(2*fftsize))

n = 1
DO k = klo,khi
  DO j = jlo,jhi
    DO i = ilo,ihi
      work(n) = n
      n = n + 1
      work(n) = n
      n = n + 1
    ENDDO
  ENDDO
ENDDO

! perform 2 FFTs

timestart = MPI_Wtime()
call fft3d_compute(fft,c_loc(work),c_loc(work),1)        ! forward FFT
CALL fft3d_compute(fft,C_LOC(work),C_LOC(work),-1)       ! inverse FFT
timestop = MPI_Wtime()

if (me == 0) then
  PRINT *,"Two",nfast,"x",nmid,"x",nslow,"FFTs on",nprocs, &
          "procs as ",npfast,"x",npmid,"x",npslow,"grid"
  PRINT *,"CPU time =",timestop-timestart,"secs"
endif

! find largest difference between initial/final values
! should be near zero

n = 1
mydiff = 0.0
DO k = klo,khi
  DO j = jlo,jhi
    DO i = ilo,ihi
      IF (abs(work(n)-n) > mydiff) mydiff = abs(work(n)-n)
      n = n + 1
      IF (abs(work(n)-n) > mydiff) mydiff = abs(work(n)-n)
      n = n + 1
    ENDDO
  ENDDO
ENDDO

CALL MPI_Allreduce(mydiff,alldiff,1,MPI_DOUBLE,MPI_MAX,world,ierr)  
IF (me == 0) PRINT *,"Max difference in initial/final values =",alldiff

! clean up

deallocate(work)
call fft3d_destroy(fft)
call MPI_Finalize()

end program simple_f90

