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
integer mode,iflag,cflag,eflag,pflag,tflag,rflag,oflag,vflag
integer seed,seedinit

integer precision
integer inxlo,inxhi,inylo,inyhi       ! initial partition of grid
integer outxlo,outxhi,outylo,outyhi   ! final partition of grid
integer nfft_in              ! # of grid pts i own in initial partition
integer nfft_out             ! # of grid pts i own in final partition
integer fftsize              ! FFT buffer size returned by FFT setup

integer tuneflag,tuneper,tuneextra
real*8  tunemax

type(c_ptr) :: fft
real(8) ::  timefft,timeinit,timesetup,timetune
real(8) :: epsmax

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

integer :: ia = 16807
integer :: im = 2147483647
real*8  :: am
integer :: iq = 127773
integer :: ir = 2836

#ifdef FFT_SINGLE
real(4), allocatable, target :: work(:)
#else
real(8), allocatable, target :: work(:)
#endif

character (len=256) :: syntax

end module data

! ---------------------------------------------------------------------
! main program
! ---------------------------------------------------------------------

program test2d_f90

use data
use iso_c_binding
use fft2d_wrap

implicit none
integer i,ierr
real*8 time1,time2

syntax = &
  "Syntax: test2d_f90 -g Nx Nx -p Px Py -n Nloop -m 0/1/2/3" // &
  c_new_line // &
  "               -i zero/step/82783 -m 0/1/2/3 -tune nper tmax extra" // &
  c_new_line // &
  "               -c point/all/combo -e pencil/brick -p array/ptr/memcpy" // &
  c_new_line // &
  "               -t -r -o -v"

#ifdef FFT_SINGLE
precision = 1;
#else
precision = 2;
#endif

am = (1.0/im)

! MPI setup

call MPI_Init(ierr)

world = MPI_COMM_WORLD
call MPI_Comm_rank(world,me,ierr)
call MPI_Comm_size(world,nprocs,ierr)

! parse command-line args

call options()

! partition FFT grid across procs, for both input and output
! create fft plan, will tune if requested
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
time2 = MPI_Wtime()
timeinit = time2 - time1

if (oflag /= 0) call output(0,"initial grid")

! perform FFTs

call MPI_Barrier(world,ierr)
time1 = MPI_Wtime()

if (mode < 2) then
  do i = 1,nloop
     call fft2d_compute(fft,c_loc(work),c_loc(work),1)
!    if (oflag /= 0) call output(1,"middle grid");
    call fft2d_compute(fft,c_loc(work),c_loc(work),-1)
  enddo
else
  do i = 1,nloop
    call fft2d_compute(fft,c_loc(work),c_loc(work),1)
  enddo   
endif

call MPI_Barrier(world,ierr)
time2 = MPI_Wtime()
timefft = time2 - time1

! validation check on result
! grid output
! timing results
! deallocate grid and plan

if (vflag == 0) call validate()
if (oflag /= 0) then
  if (mode < 2) then
    call output(0,"final grid")
  else
    call output(1,"final grid")
  endif
endif

call timing()
call deallocate_mine()
call fft2d_destroy(fft);

call MPI_Finalize(ierr)

end program test2d_f90

! ---------------------------------------------------------------------
! parse command-line options
! all options have defaults
! ---------------------------------------------------------------------

subroutine options()
use data
implicit none
integer iarg,narg
character (len=64) :: arg

! defaults

nx = 8
ny = 8
inpx = 0
inpy = 0
outpx = 0
outpy = 0
nloop = 1
iflag = zero
tuneflag = 0
mode = 0
cflag = combo
eflag = pencil
pflag = memcpy
tflag = 0
rflag = 0
oflag = 0
vflag = 0

! parse args

narg = command_argument_count()

iarg = 1
do while (iarg <= narg)
  call get_command_argument(iarg,arg)
  if (arg == '-h') then
    call error_all(syntax)
  else if (arg == '-g') then
    if (iarg+3 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') nx
    call get_command_argument(iarg+2,arg)
    read (arg,'(i10)') ny
    iarg = iarg + 3
  else if (arg == "-pin") then
    if (iarg+3 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') inpx
    call get_command_argument(iarg+2,arg)
    read (arg,'(i10)') inpy
    iarg = iarg + 3
  else if (arg == "-pout") then
    if (iarg+3 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') outpx
    call get_command_argument(iarg+2,arg) 
    read (arg,'(i10)') outpy
    iarg = iarg + 3
  else if (arg == '-n') then
    if (iarg+2 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') nloop
    iarg = iarg + 2
  else if (arg == "-i") then
    if (iarg+2 > narg+1) call error_all(syntax)
    if (arg == "zero") then
      iflag = zero
    else if (arg == "step") then
      iflag = step
    else if (arg == "index") then
      iflag = index
    else
      iflag = randominit
      ! per-processor rng seed
      call get_command_argument(iarg+1,arg)
      read (arg,'(i10)') seedinit
      seed = seedinit + me
    endif
    iarg = iarg + 2
  else if (arg == "-tune") then
    if (iarg+4 > narg+1) call error_all(syntax)
    tuneflag = 1
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') tuneper
    call get_command_argument(iarg+2,arg)
    read (arg,'(f10.3)') tunemax
    call get_command_argument(iarg+3,arg)
    read (arg,'(i10)') tuneextra
    iarg = iarg + 4
  else if (arg == "-m") then
    if (iarg+2 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    read (arg,'(i10)') mode
    iarg = iarg + 2
  else if (arg == "-c") then
    if (iarg+2 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    if (arg == "point") then
      cflag = point
    else if (arg == "all") then
      cflag = all2all
    else if (arg == "combo") then
      cflag = combo
    else 
      call error_all(syntax)
    endif
    iarg = iarg + 2
  else if (arg == "-e") then
    if (iarg+2 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    if (arg == "pencil") then
      eflag = pencil
    else if (arg == "brick") then
      eflag = brick
    else 
      call error_all(syntax)
    endif
    iarg = iarg + 2
  else if (arg == "-p") then
    if (iarg+2 > narg+1) call error_all(syntax)
    call get_command_argument(iarg+1,arg)
    if (arg == "array") then
      pflag = array
    else if (arg == "ptr") then 
      pflag = pointer
    else if (arg == "memcpy") then
      pflag = memcpy
    else 
      call error_all(syntax)
    endif
    iarg = iarg + 2
  else if (arg == "-t") then
    tflag = 1
    iarg = iarg + 1
  else if (arg == "-r") then
    rflag = 1
    iarg = iarg + 1
  else if (arg == "-o") then
    oflag = 1
    iarg = iarg + 1
  else if (arg == "-v") then
    vflag = 1
    iarg = iarg + 1
  else
    call error_all(syntax)
  endif
enddo

! sanity check on args

if (nx <= 0 .or. ny <= 0) call error_all("Invalid grid size")

if (inpx == 0 .and. inpy == 0) then
else if (inpx <= 0 .or. inpy <= 0) then
  call error_all("Invalid proc grid")
else if (inpx*inpy /= nprocs) then
  call error_all("Specified proc grid does not match nprocs")
endif

if (outpx == 0 .and. outpy == 0) then
else if (outpx <= 0 .or. outpy <= 0) then
  call error_all("Invalid proc grid")
else if (outpx*outpy /= nprocs) then
  call error_all("Specified proc grid does not match nprocs")
endif

if (nloop < 0) call error_all("Invalid nloop")
if (nloop == 0 .and. tuneflag == 0) call error_all("Invalid nloop")
if (iflag == randominit .and. seed <= 0) &
        call error_all("Invalid initialize setting")
if (mode < 0 .or. mode > 3) call error_all("Invalid FFT mode")
if (mode > 1 .and. vflag /= 0) call error_all("Cannot validate forward only FFT")

if (tuneflag /= 0 .and. tuneper <= 0) call error_all("Invalid tune nper")
if (tuneflag /= 0 .and. tunemax < 0.0) call error_all("Invalid tune tmax")
if (tuneflag /= 0 .and. (tuneextra < 0 .or. tuneextra > 1)) &
        call error_all("Invalid tune extra")
if (tuneflag /= 0 .and. rflag /= 0) call error_all("Cannot tune with remap only")

end subroutine options

! ---------------------------------------------------------------------
! partition processors across grid dimensions
! flag = in for input partitions, or out for output partitions
! if user set px,py,pz -> just return
! for in:
!   assign nprocs as bricks to 3d grid to minimize surface area per proc
!   derived from spparks domain::procs2domain_3d()
! for out:
!   assign nprocs as rectangles to xy grid to minimize surface area per proc
!   derived from spparks domain::procs2domain_2d()
! ---------------------------------------------------------------------

subroutine proc_setup(flag)
use data
implicit none
integer flag

if (flag == 0) then
  if (inpx /= 0 .or. inpy /= 0) return
  call proc2d(inpx,inpy)
endif

if (flag == 1) then
  if (outpx /= 0 .or. outpy /= 0) return
  outpx = nprocs
  outpy = 1
endif

end subroutine proc_setup

subroutine proc2d(px,py)
use data
implicit none
integer px,py
integer ipx,ipy
real*8 boxx,boxy,surf,xprd,yprd,bestsurf

xprd = nx
yprd = ny
  
bestsurf = 2.0 * (xprd+yprd)
  
! loop thru all possible factorizations of nprocs
! surf = surface area of a proc sub-domain
  
ipx = 1
do while (ipx <= nprocs)
  if (mod(nprocs,ipx) == 0) then
    ipy = nprocs/ipx
    boxx = xprd/ipx
    boxy = yprd/ipy
    surf = boxx + boxy
    if (surf < bestsurf) then
      bestsurf = surf
      px = ipx
      py = ipy
    endif
  endif
  ipx = ipx + 1
enddo
  
if (px*py /= nprocs) &
        call error_all("computed proc grid does not match nprocs")

end subroutine proc2d

! ---------------------------------------------------------------------
! partition fft grid
! once for input grid, once for output grid
! use px,py for in/out
! ---------------------------------------------------------------------

subroutine grid_setup()
use data
implicit none
integer ipx,ipy,ipz

! ipx,ipy = my position in input 3d grid of procs

ipx = mod(me,inpx)
ipy = me / inpx

! nlo,nhi = lower/upper limits of the 2d brick i own

inxlo = 1.0 * ipx * nx / inpx + 1
inxhi = 1.0 * (ipx+1) * nx / inpx

inylo = 1.0 * ipy * ny / inpy + 1
inyhi = 1.0 * (ipy+1) * ny / inpy

nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1)

! ipx,ipy,ipz = my position in output 2d grid of procs

ipx = mod(me,outpx)
ipy = me / outpx

! nlo,nhi = lower/upper limits of the 2d brick i own

outxlo = 1.0 * ipx * nx / outpx + 1
outxhi = 1.0 * (ipx+1) * nx / outpx

outylo = 1.0 * ipy * ny / outpy + 1
outyhi = 1.0 * (ipy+1) * ny / outpy

nfft_out = (outxhi-outxlo+1) * (outyhi-outylo+1)

end subroutine grid_setup

! ---------------------------------------------------------------------
! create fft plan
! ---------------------------------------------------------------------

subroutine plan()
use data
use fft2d_wrap
implicit none
integer permute,sendsize,recvsize,flag,ierr
real*8 time1,time2

call fft2d_create(world,precision,fft)
call fft2d_set(fft,"remaponly",rflag)

call fft2d_set(fft,"collective",cflag)
call fft2d_set(fft,"exchange",eflag)
call fft2d_set(fft,"pack",pflag)

if (mode == 0 .or. mode == 2) then
  permute = 0
else 
  permute = 2
endif

call MPI_Barrier(world,ierr)
time1 = MPI_Wtime()

! will use fftsize to allocate work buffer
! ignore sendsize, recvsize b/c let fft allocate remap buffers internally
! set timesetup and timetune
! reset nloop if tuning and user nloop = 0

if (tuneflag == 0) then
  call fft2d_setup(fft,nx,ny, &
          inxlo,inxhi,inylo,inyhi,outxlo,outxhi,outylo,outyhi, &
          permute,fftsize,sendsize,recvsize)
else
  flag = 0
  if (mode >= 2) flag = 1
  call fft2d_tune(fft,nx,ny, &
          inxlo,inxhi,inylo,inyhi,outxlo,outxhi,outylo,outyhi, &
          permute,fftsize,sendsize,recvsize, &
          flag,tuneper,tunemax,tuneextra)
  IF (nloop == 0) nloop = fft2d_get_int(fft,"npertrial"//c_null_char)
endif

call MPI_Barrier(world,ierr)
time2 = MPI_Wtime()

if (tuneflag == 0) then
  timesetup = time2 - time1
  timetune = 0.0
else
  timesetup = fft2d_get_double(fft,"setuptime"//c_null_char)
  timetune = time2 - time1
endif

end subroutine plan

! ---------------------------------------------------------------------
! allocate memory for fft grid
! ---------------------------------------------------------------------

subroutine allocate_mine()
use data
implicit none

allocate(work(2*fftsize))

end subroutine allocate_mine

! ---------------------------------------------------------------------
! must be called by all procs in world
! shuts down mpi and exits
! ---------------------------------------------------------------------

subroutine initialize()
use data
implicit none
integer m
integer ilocal,jlocal,iglobal,jglobal,nxlocal
real*8 random

if (iflag == zero) then
  do m = 1,2*nfft_in
    work(m) = 0.0
  enddo

else if (iflag == 1) then
  nxlocal = inxhi - inxlo + 1

  do m = 0,nfft_in-1
    ilocal = mod(m,nxlocal)
    jlocal = m / nxlocal
    iglobal = inxlo + ilocal
    jglobal = inylo + jlocal
    if (iglobal < nx/2 .and. jglobal < ny/2) then
      work(2*m+1) = 1.0
    else 
      work(2*m+1) = 0.0
      work(2*m+2) = 0.0
    endif
  enddo

else if (iflag == 2) then
  nxlocal = inxhi - inxlo + 1

  do m = 0,nfft_in-1
    ilocal = mod(m,nxlocal)
    jlocal = m / nxlocal
    iglobal = inxlo + ilocal
    jglobal = inylo + jlocal
    work(2*m+1) = jglobal + iglobal + 1
    work(2*m+2) = 0.0
  enddo

else if (iflag == 3) then
  do m = 1,2*nfft_in
    work(m) = random()
  enddo
endif

end subroutine initialize

! ---------------------------------------------------------------------
! output fft grid values
! flag = 0 for initial partition
! flag = 1 for final partition
! ---------------------------------------------------------------------

subroutine output(flag, str)
use data
implicit none
integer flag
character (len=*) :: str
integer iproc,m,tmp,ierr
integer ilocal,jlocal,iglobal,jglobal
integer nxlocal

if (me == 0) print *,str

do iproc = 0,nprocs-1
  if (me /= iproc) continue
  if (me >= 1) call MPI_Recv(tmp,0,MPI_INT,me-1,0,world,MPI_STATUS_IGNORE,ierr)

  if (flag == 0) then
    nxlocal = inxhi - inxlo + 1
    
    do m = 0,nfft_in-1
      ilocal = mod(m,nxlocal)
      jlocal = m / nxlocal
      iglobal = inxlo + ilocal
      jglobal = inylo + jlocal
      print *,"Value (",iglobal,jglobal,") on proc",me, &
              "= (",work(2*m),work(2*m+1),")"
    enddo
  else
    nxlocal = outxhi - outxlo + 1;

    do m = 0,nfft_in-1
      ilocal = mod(m,nxlocal)
      jlocal = m / nxlocal
      iglobal = outxlo + ilocal
      jglobal = outylo + jlocal
      print *,"Value (",iglobal,jglobal,") on proc",me, &
              "= (",work(2*m),work(2*m+1),")"
    enddo
  endif

  if (me < nprocs-1) call MPI_Send(tmp,0,MPI_INT,me+1,0,world,ierr);
enddo

end subroutine output

! ---------------------------------------------------------------------
! validation check for correct result
! ---------------------------------------------------------------------

subroutine validate()
use data
implicit none
integer ilocal,jlocal,iglobal,jglobal
integer nxlocal
integer m,ierr
real*8 delta,epsilon,value,newvalue
real*8 random

epsilon = 0.0

if (iflag == zero) then
  do m = 0,2*nfft_in-1
    delta = ABS(work(m+1))
    if (delta > epsilon) epsilon = delta
  enddo

else if (iflag == step) then
  nxlocal = inxhi - inxlo + 1

  do m = 0,nfft_in-1
    ilocal = mod(m,nxlocal)
    jlocal = m / nxlocal
    iglobal = inxlo + ilocal
    jglobal = inylo + jlocal
    if (iglobal < nx/2 .and. jglobal < ny/2) then
      value = 1.0
    else
      value = 0.0
    endif
    delta = ABS(work(2*m+1)-VALUE)
    if (delta > epsilon) epsilon = delta
    delta = abs(work(2*m+2))
    if (delta > epsilon) epsilon = delta
  enddo

else if (iflag == index) then
  nxlocal = inxhi - inxlo + 1

  do m = 0,nfft_in-1
    ilocal = mod(m,nxlocal)
    jlocal = m / nxlocal
    iglobal = inxlo + ilocal
    jglobal = inylo + jlocal
    value = jglobal + iglobal + 1
    delta = ABS(work(2*m+1)-VALUE)
    if (delta > epsilon) epsilon = delta
    delta = abs(work(2*m+2))
    if (delta > epsilon) epsilon = delta
  enddo

else if (iflag == randominit) then
  seed = seedinit
  do m = 0,2*nfft_in-1
    newvalue = random()
    delta = ABS(work(m+1)-newvalue)
    if (delta > epsilon) epsilon = delta
  enddo
endif

call MPI_Allreduce(epsilon,epsmax,1,MPI_DOUBLE,MPI_MAX,world,ierr)

end subroutine validate

! ---------------------------------------------------------------------
! output timing data
! ---------------------------------------------------------------------

subroutine timing()
use data
use iso_c_binding
use fft2d_wrap
implicit none
integer nfft
real (kind=8) :: onetime,nsize,log2n,floprate
integer i,nlen,ierr
real*8 time1d,time_remap;
real*8 time_remap1,time_remap2,time_remap3
real*8 time1,time2,time3,time4
integer*8 gridbytes
integer ntrial,npertrial
integer, pointer :: cflags(:) => null()
integer, pointer :: eflags(:) => null()
integer, pointer :: pflags(:) => null()
real(8), pointer :: tfft(:) => null()
real(8), pointer :: t1d(:) => null()
real(8), pointer :: tremap(:) => null()
real(8), pointer :: tremap1(:) => null()
real(8), pointer :: tremap2(:) => null()
real(8), pointer :: tremap3(:) => null()
character(c_char), pointer :: libstr(:) => null()
character(c_char), pointer :: precstr(:) => null()
type(c_ptr) :: ptr

! perform only 1d FFTs

if (tflag /= 0) then
  do i = 1,2*nfft_in
    work(i) = 0.0
  enddo

  call MPI_Barrier(world,ierr)
  time1 = MPI_Wtime()

  if (mode < 2) then
    do i = 1,nloop
      call fft2d_only_1d_ffts(fft,c_loc(work),1)
      call fft2d_only_1d_ffts(fft,c_loc(work),-1)
    enddo
  else
    do i = 1,nloop
      call fft2d_only_1d_ffts(fft,c_loc(work),1)
    enddo
  endif

  call MPI_Barrier(world,ierr)
  time2 = MPI_Wtime()
  time1d = time2 - time1
endif

! perform all remaps

if (tflag /= 0) then
  do i = 1,2*nfft_in
    work(i) = 0.0
  enddo

  call MPI_Barrier(world,ierr)
  time1 = MPI_Wtime()

  if (mode < 2) then
    do i = 1,nloop
      call fft2d_only_remaps(fft,c_loc(work),c_loc(work),1)
      call fft2d_only_remaps(fft,c_loc(work),c_loc(work),-1)
    enddo
  else
    do i = 1,nloop
      call fft2d_only_remaps(fft,c_loc(work),c_loc(work),1)
    enddo
  endif

  call MPI_Barrier(world,ierr)
  time2 = MPI_Wtime()
  time_remap = time2 - time1
endif

! perform only single remaps

if (tflag /= 0) then
  do i = 1,2*nfft_in
    work(i) = 0.0
  enddo

  call MPI_Barrier(world,ierr)
  time1 = MPI_Wtime()

  if (mode < 2) then
    do i = 1,nloop
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),1,1)
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),-1,1)
    enddo
  else
    do i = 1,nloop
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),1,1)
    enddo
  endif

  call MPI_Barrier(world,ierr)
  time2 = MPI_Wtime()
  time_remap1 = time2 - time1

  if (mode < 2) then
    do i = 1,nloop
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),1,2)
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),-1,2)
    enddo
  else
    do i = 1,nloop
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),1,2)
    enddo
  endif

  call MPI_Barrier(world,ierr)
  time3 = MPI_Wtime()
  time_remap2 = time3 - time2

  if (mode < 2) then
    do i = 1,nloop
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),1,3)
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),-1,3)
    enddo
  else
    do i = 1,nloop
      call fft2d_only_one_remap(fft,c_loc(work),c_loc(work),1,3)
    enddo
  endif

  call MPI_Barrier(world,ierr)
  time4 = MPI_Wtime()
  time_remap3 = time4 - time3
endif

! stats output
! nfft = 2x larger for modes 0,1

if (mode < 2) then
  nfft = 2*nloop
else 
  nfft = nloop
endif

onetime = timefft/nfft
nsize = nx * ny
log2n = log(nsize)/log(2.0)
floprate = 5.0 * nsize * log2n / onetime / (1024*1024*1024)

#ifdef FFT_SINGLE
gridbytes = 4 * 2*fftsize
#else
gridbytes = 8 * 2*fftsize
#endif

if (me == 0) then
  ptr = fft2d_get_string(fft,"fft1d"//c_null_char,nlen)
  call c_f_pointer(ptr,libstr,[nlen])
  ptr = fft2d_get_string(fft,"precision"//c_null_char,nlen)
  call c_f_pointer(ptr,precstr,[nlen])

  print *,"2d FFTs with ",libstr," library, precision = ",precstr
  print *,"Grid size:",nx,ny
  print *,"  initial proc grid:",inpx,inpy
  print *,"  x pencil proc grid:", &
          fft2d_get_int(fft,"npfast1"//c_null_char), &
          fft2d_get_int(fft,"npfast2"//c_null_char)
  print *,"  y pencil proc grid:", &
          fft2d_get_int(fft,"npslow1"//c_null_char), &
          fft2d_get_int(fft,"npslow2"//c_null_char)
  print *,"  2d brick proc grid:", &
          fft2d_get_int(fft,"npbrick1"//c_null_char), &
          fft2d_get_int(fft,"npbrick2"//c_null_char)
  print *,"  final proc grid:",outpx,outpy

  if (tuneflag /= 0) then
    ntrial = fft2d_get_int(fft,"ntrial"//c_null_char)
    npertrial = fft2d_get_int(fft,"npertrial"//c_null_char)
    print *,"Tuning trials & iterations:",ntrial,npertrial
    ptr = fft2d_get_int_vector(fft,"cflags"//c_null_char,nlen)
    call c_f_pointer(ptr,cflags,[nlen])
    ptr = fft2d_get_int_vector(fft,"eflags"//c_null_char,nlen)
    call c_f_pointer(ptr,eflags,[nlen])
    ptr = fft2d_get_int_vector(fft,"pflags"//c_null_char,nlen)
    call c_f_pointer(ptr,pflags,[nlen])
    ptr = fft2d_get_double_vector(fft,"tfft"//c_null_char,nlen)
    call c_f_pointer(ptr,tfft,[nlen])
    ptr = fft2d_get_double_vector(fft,"t1d"//c_null_char,nlen)
    call c_f_pointer(ptr,t1d,[nlen])
    ptr = fft2d_get_double_vector(fft,"tremap"//c_null_char,nlen)
    call c_f_pointer(ptr,tremap,[nlen])
    ptr = fft2d_get_double_vector(fft,"tremap1"//c_null_char,nlen)
    call c_f_pointer(ptr,tremap1,[nlen])
    ptr = fft2d_get_double_vector(fft,"tremap2"//c_null_char,nlen)
    call c_f_pointer(ptr,tremap2,[nlen])
    ptr = fft2d_get_double_vector(fft,"tremap3"//c_null_char,nlen)
    call c_f_pointer(ptr,tremap3,[nlen])
    do i = 1,ntrial
      print *,"  coll exch pack 2dFFT 1dFFT remap r1 r2 r3:", &
              cflags(i),eflags(i),pflags(i),tfft(i),t1d(i),tremap(i), &
              tremap1(i),tremap2(i),tremap3(i)
    enddo
  endif

  if (mode == 0) then
    print *,nloop,"forward and",nloop,"back FFTs on",nprocs,"procs"
  else if (mode == 1) then
    print *,nloop,"forward and",nloop,"back convolution FFTs on",nprocs,"procs"
  else if (mode == 2) then
    print *,nloop,"forward FFTs on",nprocs,"procs"
  else if (mode == 3) then
    print *,nloop,"forward convolution FFTs on",nprocs,"procs"
  endif

  print *,"Collective, exchange, pack methods:", &
          fft2d_get_int(fft,"collective"//c_null_char), &
          fft2d_get_int(fft,"exchange"//c_null_char), &
          fft2d_get_int(fft,"pack"//c_null_char)
  print *,"Memory usage (per-proc) for FFT grid =", &
          1.0*gridbytes / 1024/1024,"MBytes"
  print *,"Memory usage (per-proc) by fftMPI =", &
          1.0*fft2d_get_int64(fft,"memusage"//c_null_char) / 1024/1024,"MBytes"
  
  if (vflag /= 0) print *,"Max error =",epsmax
  if (tuneflag /= 0) then
    print *,"Initialize grid =",timeinit-timesetup,"secs"
  else 
    print *,"Initialize grid =",timeinit-timetune,"secs"
  endif
  print *,"FFT setup =",timesetup,"secs"
  print *,"FFT tune =",timetune,"secs"
  print *,"Time for 2d FFTs =",timefft,"secs"
  print *,"  time/fft2d =",onetime,"secs"
  print *,"  flop rate for 2d FFTs =",floprate,"Gflops"
  if (tflag /= 0) then
    print *,"Time for 1d FFTs only =",time1d,"secs"
    print *,"  time/fft1d =",time1d/nfft,"secs"
    print *,"  fraction of time in 1d FFTs =",time1d/timefft
  endif
  if (tflag /= 0) then
    print *,"Time for remaps only =",time_remap,"secs"
    print *,"  fraction of time in remaps =",time_remap/timefft
    print *,"Time for remap #1 =",time_remap1,"secs"
    print *,"  fraction of time in remap #1 =",time_remap1/timefft
    print *,"Time for remap #2 =",time_remap2,"secs"
    print *,"  fraction of time in remap #2 =",time_remap2/timefft
    print *,"Time for remap #3 =",time_remap3,"secs"
    print *,"  fraction of time in remap #3 =",time_remap3/timefft
  endif
endif

end subroutine timing

! ---------------------------------------------------------------------
! deallocate memory for fft grid
! ---------------------------------------------------------------------

subroutine deallocate_mine()
use data
implicit none

deallocate(work)

end subroutine deallocate_mine

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! utility functions
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! must be called by all procs in world
! shuts down MPI and exits
! ---------------------------------------------------------------------

subroutine error_all(str)
use data
implicit none
character (len=*) :: str
integer ierr

call MPI_Barrier(world,ierr)
if (me == 0) print *,"ERROR: ",str
call MPI_Finalize(ierr)
call exit()

end subroutine error_all

! ----------------------------------------------------------------------
! simple park rng
! pass in non-zero seed
! ----------------------------------------------------------------------

function random()
use data
implicit none
integer k
real*8 ans,random

k = seed/iq
seed = ia*(seed-k*iq) - ir*k
if (seed < 0) seed = seed + im
ans = am*seed
random = ans
return

end function random
