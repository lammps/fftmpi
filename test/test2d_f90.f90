! test driver on 2d FFT from FFT library
! for benchmarking, timing purposes
! see test2d.cpp for command-line args

! ---------------------------------------------------------------------
! common data
! ---------------------------------------------------------------------

module data
use iso_c_binding
implicit none

include 'mpif.h'

integer world
integer me,nprocs

integer nx,ny
integer inpx,inpy,outpx,outpy
integer nloop
INTEGER mode,iflag,cflag,eflag,pflag,tflag,rflag,oflag,vflag
integer seed

integer precision
integer inxlo,inxhi,inylo,inyhi       ! initial partition of grid
integer outxlo,outxhi,outylo,outyhi   ! final partition of grid
integer nfft_in              ! # of grid pts I own in initial partition
integer nfft_out             ! # of grid pts I own in final partition
integer fftsize              ! FFT buffer size returned by FFT setup

integer tuneflag,tuneper,tuneextra
REAL*8  tunemax

TYPE(C_ptr) :: fft
REAL*8 :: time2d,timeinit

integer :: zero = 0
integer :: step = 1
integer :: index = 2
integer :: randominit = 3
integer :: point = 0
integer :: all2all = 1
integer :: combo = 2
integer :: pencil = 0
integer :: brick = 1
integer :: array = 0
integer :: pointer = 1
integer :: memcpy = 2
integer :: in = 0
integer :: out = 1

integer :: IA = 16807
integer :: IM = 2147483647
REAL*8  :: AM
integer :: IQ = 127773
integer :: IR = 2836

#ifdef FFT_SINGLE
REAL(4), ALLOCATABLE, target :: work(:)
#else
REAL(8), ALLOCATABLE, target :: work(:)
#endif

end module data

! ---------------------------------------------------------------------
! main program
! ---------------------------------------------------------------------

program test2d_f90

use data
use iso_c_binding
use fft2d_wrap

implicit none
INTEGER i,ierr
REAL*8 time1,time2

#ifdef FFT_SINGLE
precision = 1;
#else
precision = 2;
#endif

AM = (1.0/IM)

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
ALLOCATE(work(2*fftsize))
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
    call fft2d_compute(fft,c_loc(work),c_loc(work),1)
!    if (oflag /= 0) call output(1,"Middle grid");
    call fft2d_compute(fft,c_loc(work),c_loc(work),-1)
  enddo
else
  do i = 1,nloop
    call fft2d_compute(fft,c_loc(work),c_loc(work),1)
  enddo   
endif

call MPI_Barrier(world,ierr)
time2 = mpi_wtime()
time2d = time2 - time1

! validation check on result
! grid output
! timing results
! deallocate grid and plan

if (vflag == 0) call validate()
IF (oflag /= 0) THEN
  if (mode < 2) then
    call output(0,"Final grid")
  else
    call output(1,"Final grid")
  endif
endif

call timing()
deallocate(work)
call fft2d_destroy(fft);

call MPI_Finalize(ierr)

end program test2d_f90

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

syntax = "Syntax: test2d_f90 -g Nx Ny -n Niter"

! defaults

nx = 8
ny = 8
inpx = 0
inpy = 0
outpx = 0
outpy = 0
nloop = 1
iflag = ZERO
tuneflag = 0
mode = 0
cflag = COMBO
eflag = PENCIL
pflag = MEMCPY
tflag = 0
rflag = 0
oflag = 0
vflag = 0

! parse args

narg = command_argument_count()

iarg = 1
do while (iarg <= narg)
  call get_command_argument(iarg,arg)
  IF (arg == '-h') THEN
    call error_all(syntax)
  else if (arg == '-g') then
    IF (iarg+2 > narg) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') nx
    call get_command_argument(iarg+2,arg)
    read (arg,'(i10)') ny
    iarg = iarg + 3
  ELSE IF (arg == "-pin") then
    IF (iarg+3 > narg) error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') inpx
    call get_command_argument(iarg+2,arg)
    read (arg,'(i10)') inpy
    iarg += 3
  else if (args == "-pout") then
    IF (iarg+3 > narg) error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') outpx
    call get_command_argument(iarg+2,arg) 
    read (arg,'(i10)') outpy
    iarg += 3
  else if (arg == '-n') then
    IF (iarg+1 > narg) CALL error_all(syntax)
    CALL GET_COMMAND_ARGUMENT(iarg+1,arg)
    READ (arg,'(i10)') nloop
    iarg = iarg + 2
  ELSE IF (arg == "-i") THEN
    IF (iarg+2 > narg) error_all(syntax)
    IF (arg == "zero") THEN
      flag = ZERO
    ELSE IF (arg == "step") then
      flag = STEP
    ELSE IF (arg == "index") then
      flag = INDEX
    ELSE
      iflag = RANDOM
      ! per-processor RNG seed
      seed = seedinit = atoi(args[iarg+1]) + me
    ENDIF
    iarg += 2
  ELSE IF (arg == "-tune") THEN
    IF (iarg+4 > narg) error_all(syntax)
    tuneflag = 1
    tuneper = atoi(args[iarg+1])
    tunemax = atof(args[iarg+2])
    tuneextra = atoi(args[iarg+3])
    iarg += 4
  ELSE IF (arg == "-m") THEN
    IF (iarg+2 > narg) error_all(syntax)
    mode = atoi(args[iarg+1])
    iarg += 2
  ELSE IF (arg == "-c") THEN
    IF (iarg+2 > narg) error_all(syntax)
    IF (arg == "point") THEN
      flag = POINT
    ELSE IF (arg == "all") then
      flag = ALL2ALL
    ELSE IF (arg == "combo") then
      flag = COMBO
    ELSE error_all(syntax)
    ENDIF
    iarg += 2
  ELSE IF (arg == "-e") THEN
    IF (iarg+2 > narg) error_all(syntax)
    IF (arg == "pencil") then
      flag = PENCIL
    ELSE IF (arg == "brick") then
      flag = BRICK
    ELSE error_all(syntax)
    endif
    iarg += 2
  ELSE IF (arg == "-p") THEN
    IF (iarg+2 > narg) error_all(syntax)
    IF (arg == "array") THEN
      flag = ARRAY
    ELSE IF (arg == "ptr") THEN 
      flag = POINTER
    ELSE IF (arg == "memcpy") THEN
      flag = MEMCPY
    ELSE error_all(syntax)
    ENDIF
    iarg += 2
  ELSE IF (arg == "-t") THEN
    tflag = 1
    iarg += 1
  ELSE IF (arg == "-r") THEN
    rflag = 1
    iarg += 1
  ELSE IF (arg == "-o") THEN
    oflag = 1
    iarg += 1
  ELSE IF (arg == "-v") THEN
    vflag = 1
    iarg += 1
  else
    call error_all(syntax)
  endif
enddo

! sanity check on args

if (nx <= 0 .or. ny <= 0) call error_all("Invalid grid size")
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
  IF (inpx /= 0 .OR. inpy /= 0) RETURN
  call proc2d(inpx,inpy)
endif

IF (flag == 1) then
  IF (outpx /= 0 .OR. outpy /= 0) RETURN
  outpx = nprocs
  outpy = 1
endif

end subroutine proc_setup

SUBROUTINE proc2d(px,py)
use data
implicit none
INTEGER px,py
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
  
IF (px*py /= nprocs) &
        CALL error_all("Computed proc grid does not match nprocs")

end subroutine proc2d

! ---------------------------------------------------------------------
! partition FFT grid
! once for input grid, once for output grid
! use Px,Py for in/out
! ---------------------------------------------------------------------

SUBROUTINE grid_setup()
use data
implicit none
INTEGER ipx,ipy,ipz

! ipx,ipy = my position in input 3d grid of procs

ipx = MOD(me,inpx)
ipy = me / inpx

! nlo,nhi = lower/upper limits of the 2d brick I own

inxlo = 1.0 * ipx * nx / inpx + 1
inxhi = 1.0 * (ipx+1) * nx / inpx

inylo = 1.0 * ipy * ny / inpy + 1
inyhi = 1.0 * (ipy+1) * ny / inpy

nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1)

! ipx,ipy,ipz = my position in output 2d grid of procs

ipx = MOD(me,outpx)
ipy = me / outpx

! nlo,nhi = lower/upper limits of the 2d brick I own

outxlo = 1.0 * ipx * nx / outpx + 1
outxhi = 1.0 * (ipx+1) * nx / outpx

outylo = 1.0 * ipy * ny / outpy + 1
outyhi = 1.0 * (ipy+1) * ny / outpy

nfft_out = (outxhi-outxlo+1) * (outyhi-outylo+1)

end subroutine grid_setup

! ---------------------------------------------------------------------
! must be called by all procs in world
! shuts down MPI and exits
! ---------------------------------------------------------------------

SUBROUTINE plan()
use data
use fft2d_wrap
implicit none
INTEGER permute
integer sendsize,recvsize

call fft2d_create(world,precision,fft)

permute = 0
call fft2d_setup(fft,nx,ny, &
    inxlo,inxhi,inylo,inyhi, &
    outxlo,outxhi,outylo,outyhi, &
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
INTEGER ilocal,jlocal,iglobal,jglobal,nxlocal
real*8 random

if (iflag == 0) then
  DO m = 1,2*nfft_in
    work(m) = 0.0
  enddo

else if (iflag == 1) then
  nxlocal = inxhi - inxlo + 1
  DO m = 0,nfft_in-1
    ilocal = MOD(m,nxlocal)
    jlocal = m / nxlocal
    iglobal = inxlo + ilocal
    jglobal = inylo + jlocal
    IF (iglobal < nx/2 .and. jglobal < ny/2) THEN
      work(2*m+1) = 1.0
    ELSE 
      work(2*m+1) = 0.0
      work(2*m+2) = 0.0
    ENDIF
  ENDDO

else if (iflag == 2) then
  nxlocal = inxhi - inxlo + 1
  DO m = 0,nfft_in-1
    ilocal = MOD(m,nxlocal)
    jlocal = m / nxlocal
    iglobal = inxlo + ilocal
    jglobal = inylo + jlocal
    work(2*m+1) = jglobal + iglobal + 1
    work(2*m+2) = 0.0
  ENDDO

else if (iflag == 3) then
  DO m = 1,2*nfft_in
    work(m) = random()
  ENDDO
endif

end subroutine initialize

! ----------------------------------------------------------------------
! output FFT grid values
! flag = 0 for initial partition
! flag = 1 for final partition
! ----------------------------------------------------------------------

subroutine output(flag, str)
use data
implicit none
integer flag
CHARACTER(*) str

end subroutine output

! ----------------------------------------------------------------------
! validation check for correct result
! ----------------------------------------------------------------------

subroutine validate()
use data
implicit none

end subroutine validate

! ---------------------------------------------------------------------
! output timing data
! ---------------------------------------------------------------------

SUBROUTINE timing()
use data
implicit none
integer nfft
REAL (kind=8) :: onetime,nsize,log2n,floprate

nfft = 2*nloop
onetime = time2d/nfft
nsize = 1.0 * nx * ny
log2n = log(nsize)/log(2.0)
floprate = 5.0 * nsize * log2n / onetime / (1024*1024*1024)

if (me == 0) then
!  printf("2d FFTs with %s library, precision = %s\n",
!    fft->fft1d,fft->PRECISION)

  PRINT *,"Grid size:",nx,ny

  PRINT *,nloop,"forward and ",nloop,"back FFTs on",nprocs,"procs"
  PRINT *,"Setup time =",timeinit,"secs"
  PRINT *,"Time for 2d FFTs =",time2d,"secs"
  PRINT *,"  time/fft2d =",onetime,"secs"
  PRINT *,"  flop rate for 2d FFTs =",floprate,"Gflops"
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
if (me == 0) print *,"ERROR:",trim(str)
CALL MPI_Finalize(ierr)
call exit()

end subroutine error_all

! ----------------------------------------------------------------------
! simple Park RNG
! pass in non-zero seed
! ----------------------------------------------------------------------

function random()
use data
implicit none
integer k
REAL*8 ans,random

k = seed/IQ
seed = IA*(seed-k*IQ) - IR*k
IF (seed < 0) seed = seed + IM
ans = AM*seed
random = ans
return

end function random
