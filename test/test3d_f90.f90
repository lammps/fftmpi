! test driver on 3d FFT from FFT library
! for benchmarking, timing purposes
! see test3d.cpp for command-line args

! ---------------------------------------------------------------------
! common data
! ---------------------------------------------------------------------

module data
use iso_c_binding
implicit none

include 'mpif.h'

integer world
integer me,nprocs

integer nx,ny,nz
integer inpx,inpy,inpz,outpx,outpy,outpz
INTEGER nloop
INTEGER mode,oflag

integer precision
integer inxlo,inxhi,inylo,inyhi,inzlo,inzhi  ! initial partition of grid
integer outxlo,outxhi,outylo,outyhi,outzlo,outzhi   ! final partition of grid
integer nfft_in              ! # of grid pts I own in initial partition
integer nfft_out             ! # of grid pts I own in final partition
integer fftsize              ! FFT buffer size returned by FFT setup

TYPE(C_ptr) :: fft
REAL(8) :: time3d,timeinit

#ifdef FFT_SINGLE
REAL(4), ALLOCATABLE, target :: work(:)
#else
REAL(8), ALLOCATABLE, target :: work(:)
#endif

end module data

! ---------------------------------------------------------------------
! main program
! ---------------------------------------------------------------------

program test3d_f90

use data
use iso_c_binding
use fft3d_wrap

implicit none
INTEGER i,ierr
REAL*8 time1,time2

#ifdef FFT_SINGLE
precision = 1;
#else
precision = 2;
#endif

! MPI setup

call MPI_Init(ierr)

world = MPI_COMM_WORLD
call MPI_Comm_rank(world,me,ierr)
call MPI_Comm_size(world,nprocs,ierr)

! parse command-line args

call options()

! partition FFT grid across procs, for both input and output
! create FFT plan, will tune if requested
! allocate grid
! initialize FFT grid
! grid output

call MPI_Barrier(world,ierr)
time1 = MPI_Wtime()

call proc_setup(0)
call proc_setup(1)
call grid_setup()
call plan()
call allocate_mine()
call initialize()

call MPI_Barrier(world,ierr)
time2 = mpi_wtime()
timeinit = time2 - time1

IF (oflag /= 0) CALL output(0,"Initial grid")

! perform FFTs

call MPI_Barrier(world,ierr)
time1 = MPI_Wtime()

if (mode < 2) then
  do i = 1,nloop
    call fft3d_compute(fft,c_loc(work),c_loc(work),1)
!    if (oflag /= 0) call output(1,"Middle grid")
    call fft3d_compute(fft,c_loc(work),c_loc(work),-1)
  enddo
else
  do i = 1,nloop
    call fft3d_compute(fft,c_loc(work),c_loc(work),1)
  enddo
endif

call MPI_Barrier(world,ierr)
time2 = mpi_wtime()
time3d = time2 - time1

! validation check on result
! grid output
! timing results
! deallocate grid and plan

call validate()
IF (oflag /= 0) THEN
   if (mode < 2) then
      call output(0,"Final grid")
   else
      call output(1,"Final grid")
   endif
endif
  
call timing()
call deallocate_mine()
call fft3d_destroy(fft)

call MPI_Finalize(ierr)

end program test3d_f90

! ---------------------------------------------------------------------
! parse command-line options
! all options have defaults
! ---------------------------------------------------------------------

SUBROUTINE options()
use data
implicit none
INTEGER iarg,narg
character (len=64) :: arg
character (len=128) :: syntax

syntax = "Syntax: test3d_f90 -g Nx Ny Nz -n Niter"

! defaults

nx = 8
ny = 8
nz = 8
inpx = 0
inpy = 0
inpz = 0
outpx = 0
outpy = 0
outpz = 0
nloop = 1

! parse args

narg = command_argument_count()

iarg = 1
do while (iarg < narg)
  call get_command_argument(iarg,arg)
  if (arg == '-h') then
    call error_all(syntax)
  else if (arg == '-n') then
    IF (iarg+1 > narg) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') nloop
    iarg = iarg + 2
  else if (arg == '-g') then
    IF (iarg+3 > narg) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') nx
    call get_command_argument(iarg+2,arg)
    read (arg,'(i10)') ny
    call get_command_argument(iarg+3,arg)
    read (arg,'(i10)') nz
    iarg = iarg + 4
  else
    call error_all(syntax)
  endif
enddo

! sanity check on args

if (nx <= 0 .or. ny <= 0 .or. nz <= 0) call error_all("Invalid grid size")
if (nloop <= 0) call error_all("Invalid Niter")

end subroutine options

! ---------------------------------------------------------------------
! partition processors across grid dimensions
! flag = IN for input partitions, or OUT for output partitions
! if user set Px,Py,Pz -> just return
! for IN:
!   assign nprocs as bricks to 3d grid to minimize surface area per proc
!   derived from SPPARKS Domain::procs2domain_3d()
! for OUT:
!   assign nprocs as rectangles to xy grid to minimize surface area per proc
!   derived from SPPARKS Domain::procs2domain_2d()
! ---------------------------------------------------------------------

subroutine proc_setup(flag)
use data
implicit none
INTEGER flag

IF (flag == 0) then
  IF (inpx /= 0 .OR. inpy /= 0 .OR. inpz /= 0) RETURN
  call proc3d(inpx,inpy,inpz)
endif

IF (flag == 1) then
  IF (outpx /= 0 .OR. outpy /= 0 .OR. outpz /= 0) RETURN
  call proc3d(outpx,outpy,outpz)
endif

end subroutine proc_setup

SUBROUTINE proc3d(px,py,pz)
use data
implicit none
INTEGER px,py,pz
integer ipx,ipy,ipz,nremain
REAL*8 boxx,boxy,boxz,surf
REAL*8 xprd,yprd,zprd,bestsurf

xprd = nx
yprd = ny
zprd = nz
  
bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd)
  
! loop thru all possible factorizations of nprocs
! surf = surface area of a proc sub-domain
  
ipx = 1
DO WHILE (ipx <= nprocs)
  IF (MOD(nprocs,ipx) == 0) THEN
    nremain = nprocs/ipx
    ipy = 1
    DO WHILE (ipy <= nremain)
      IF (MOD(nremain,ipy) == 0) THEN
        ipz = nremain/ipy
        boxx = xprd/ipx
        boxy = yprd/ipy
        boxz = zprd/ipz
        surf = boxx*boxy + boxy*boxz + boxz*boxx
        IF (surf < bestsurf) THEN
          bestsurf = surf
          px = ipx
          py = ipy
          pz = ipz
        ENDIF
      ENDIF
      ipy = ipy + 1
    ENDDO
  ENDIF
  ipx = ipx + 1
ENDDO
  
IF (px*py*pz /= nprocs) &
        CALL error_all("Computed proc grid does not match nprocs")

end subroutine proc3d

SUBROUTINE proc2d(px,py,pz)
use data
implicit none
INTEGER px,py,pz
integer ipx,ipy
REAL*8 boxx,boxy,surf,xprd,yprd,bestsurf

xprd = nx
yprd = ny
  
bestsurf = 2.0 * (xprd+yprd)
  
! loop thru all possible factorizations of nprocs
! surf = surface area of a proc sub-domain
  
ipx = 1
do while (ipx <= nprocs)
  IF (MOD(nprocs,ipx) == 0) then
    ipy = nprocs/ipx
    boxx = xprd/ipx
    boxy = yprd/ipy
    surf = boxx + boxy
    IF (surf < bestsurf) then
      bestsurf = surf
      px = ipx
      py = ipy
    endif
  endif
  ipx = ipx + 1
enddo
  
pz = 1
IF (px*py*pz /= nprocs) &
        CALL error_all("Computed proc grid does not match nprocs")

end subroutine proc2d

! ---------------------------------------------------------------------
! partition FFT grid
! once for input grid, once for output grid
! use Px,Py,Pz for in/out
! ---------------------------------------------------------------------

SUBROUTINE grid_setup()
use data
implicit none
INTEGER ipx,ipy,ipz

! ipx,ipy,ipz = my position in input 3d grid of procs

ipx = MOD(me,inpx)
ipy = MOD(me/inpx,inpy)
ipz = me / (inpx*inpy)

! nlo,nhi = lower/upper limits of the 3d brick I own

inxlo = 1.0 * ipx * nx / inpx + 1
inxhi = 1.0 * (ipx+1) * nx / inpx

inylo = 1.0 * ipy * ny / inpy + 1
inyhi = 1.0 * (ipy+1) * ny / inpy

inzlo = 1.0 * ipz * nz / inpz + 1
inzhi = 1.0 * (ipz+1) * nz / inpz

nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1)

! ipx,ipy,ipz = my position in output 3d grid of procs

ipx = MOD(me,outpx)
ipy = MOD(me/outpx,outpy)
ipz = me / (outpx*outpy)

! nlo,nhi = lower/upper limits of the 3d brick I own

outxlo = 1.0 * ipx * nx / outpx + 1
outxhi = 1.0 * (ipx+1) * nx / outpx

outylo = 1.0 * ipy * ny / outpy + 1
outyhi = 1.0 * (ipy+1) * ny / outpy

outzlo = 1.0 * ipz * nz / outpz + 1
outzhi = 1.0 * (ipz+1) * nz / outpz

nfft_out = (outxhi-outxlo+1) * (outyhi-outylo+1) * (outzhi-outzlo+1)

end subroutine grid_setup

! ---------------------------------------------------------------------
! must be called by all procs in world
! shuts down MPI and exits
! ---------------------------------------------------------------------

SUBROUTINE plan()
use data
use fft3d_wrap
implicit none
INTEGER permute,sendsize,recvsize

call fft3d_create(world,precision,fft)

permute = 0
call fft3d_setup(fft,nx,ny,nz, &
    inxlo,inxhi,inylo,inyhi,inzlo,inzhi, &
    outxlo,outxhi,outylo,outyhi,outzlo,outzhi, &
    permute,fftsize,sendsize,recvsize)

end subroutine plan

! ---------------------------------------------------------------------
! must be called by all procs in world
! shuts down MPI and exits
! ---------------------------------------------------------------------

SUBROUTINE initialize()
use data
implicit none
integer m

DO m = 1,2*nfft_in
  work(m) = 0.0
enddo

end subroutine initialize

! ---------------------------------------------------------------------
! output timing data
! ---------------------------------------------------------------------

SUBROUTINE timing()
use data
implicit none
integer nfft
REAL (kind=8) :: onetime,nsize,log2n,floprate

nfft = 2*nloop
onetime = time3d/nfft
nsize = 1.0 * nx * ny * nz
log2n = log(nsize)/log(2.0)
floprate = 5.0 * nsize * log2n / onetime / (1024*1024*1024)

if (me == 0) then
  PRINT *,"Grid size:",nx,ny,nz
  PRINT *,nloop,"forward and",nloop,"back FFTs on",nprocs,"procs"
  PRINT *,"Setup time =",timeinit,"secs"
  PRINT *,"Time for 3d FFTs =",time3d,"secs"
  PRINT *,"  time/fft3d =",onetime,"secs"
  PRINT *,"  flop rate for 3d FFTs =",floprate,"Gflops"
endif

end subroutine timing

! ---------------------------------------------------------------------
! must be called by all procs in world
! shuts down MPI and exits
! ---------------------------------------------------------------------

SUBROUTINE error_all(str)
use data
implicit none
CHARACTER (len=*) :: str
INTEGER ierr

CALL MPI_Barrier(world,ierr)
if (me == 0) print *,"ERROR:",str
CALL MPI_Finalize(ierr)
call exit()

end subroutine error_all
