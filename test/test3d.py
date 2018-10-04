#!/usr/bin/env python

# test driver on 3d FFT from FFT library
# for benchmarking, timing purposes
# see test2d.cpp for command-line args

import sys,math
import numpy as np
from mpi4py import MPI
from fftmpi import FFT3dMPI

IA = 16807
IM = 2147483647
AM = 1.0/IM
IQ = 127773
IR = 2836

syntax = \
  "Syntax: test3d.py -g Nx Ny Nz -p Px Py Pz -n Nloop -m 0/1/2/3\n" + \
  "               -i zero/step/82783 -m 0/1/2/3 -tune nper tmax extra\n" + \
  "               -c point/all/combo -e pencil/brick -p array/ptr/memcpy\n" + \
  "               -t -r -o -v"

ZERO,STEP,INDEX,RANDOM = range(4)
POINT,ALL2ALL,COMBO = range(3)
PENCIL,BRICK = range(2)
ARRAY,POINTER,MEMCPY = range(3)
IN,OUT = range(2)

# NOTE: how to set this at run-time

precision = 2

# ----------------------------------------------------------------------
# must be called by all procs in world
# ----------------------------------------------------------------------

def error_all(txt):
  if world.rank == 0: print "ERROR:",txt
  sys.exit()

# ----------------------------------------------------------------------
# simple Park RNG
# pass in non-zero seed
# ----------------------------------------------------------------------

def random():
  global seed
  k = seed/IQ
  seed = IA*(seed-k*IQ) - IR*k
  if seed < 0: seed += IM
  ans = AM*seed
  return ans


# ----------------------------------------------------------------------
# parse command-line options
# all options have defaults
# ----------------------------------------------------------------------

def options():
  global nx,ny,nz,inpx,inpy,inpz,outpx,outpy,outpz,nloop,iflag,seed,seedinit
  global tuneflag,tuneper,tunemax,tuneextra,mode
  global cflag,eflag,pflag,tflag,rflag,oflag,vflag
  
  args = sys.argv
  narg = len(args)

  # defaults

  nx = ny = nz = 8
  inpx = inpy = inpz = 0
  outpx = outpy = outpz = 0
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

  # parse args

  iarg = 1
  while iarg < narg:
    if args[iarg] == "-h": error_all(syntax)
    elif args[iarg] == "-g":
      if iarg+4 > narg: error_all(syntax)
      nx = int(args[iarg+1])
      ny = int(args[iarg+2])
      nz = int(args[iarg+3])
      iarg += 4
    elif args[iarg] == "-pin":
      if iarg+4 > narg: error_all(syntax)
      inpx = int(args[iarg+1])
      inpy = int(args[iarg+2])
      inpz = int(args[iarg+3])
      iarg += 4
    elif args[iarg] == "-pout":
      if iarg+4 > narg: error_all(syntax)
      outpx = int(args[iarg+1])
      outpy = int(args[iarg+2])
      outpz = int(args[iarg+3])
      iarg += 4
    elif args[iarg] == "-n":
      if iarg+2 > narg: error_all(syntax)
      nloop = int(args[iarg+1])
      iarg += 2
    elif args[iarg] == "-i":
      if iarg+2 > narg: error_all(syntax)
      if args[iarg+1] == "zero": iflag = ZERO
      elif args[iarg+1] == "step": iflag = STEP
      elif args[iarg+1] == "index": iflag = INDEX
      else:
        iflag = RANDOM
        # per-processor RNG seed
        seed = seedinit = int(args[iarg+1]) + me
      iarg += 2
    elif args[iarg] == "-tune":
      if iarg+4 > narg: error_all(syntax)
      tuneflag = 1
      tuneper = int(args[iarg+1])
      tunemax = float(args[iarg+2])
      tuneextra = int(args[iarg+3])
      iarg += 4
    elif args[iarg] == "-m":
      if iarg+2 > narg: error_all(syntax)
      mode = int(args[iarg+1])
      iarg += 2
    elif args[iarg] == "-c":
      if iarg+2 > narg: error_all(syntax)
      if args[iarg+1] == "point": cflag = POINT
      elif args[iarg+1] == "all": cflag = ALL2ALL
      elif args[iarg+1] == "combo": cflag = COMBO
      else: error_all(syntax)
      iarg += 2
    elif args[iarg] == "-e":
      if iarg+2 > narg: error_all(syntax)
      if args[iarg+1] == "pencil": eflag = PENCIL
      elif args[iarg+1] == "brick": eflag = BRICK
      else: error_all(syntax)
      iarg += 2
    elif args[iarg] == "-p":
      if iarg+2 > narg: error_all(syntax)
      if args[iarg+1] == "array": pflag = ARRAY
      elif args[iarg+1] == "ptr": pflag = POINTER
      elif args[iarg+1] == "memcpy": pflag = MEMCPY
      else: error_all(syntax)
      iarg += 2
    elif args[iarg] == "-t":
      tflag = 1
      iarg += 1
    elif args[iarg] == "-r":
      rflag = 1
      iarg += 1
    elif args[iarg] == "-o":
      oflag = 1
      iarg += 1
    elif args[iarg] == "-v":
      vflag = 1
      iarg += 1
    else: error_all(syntax)

  # sanity check on args

  if nx <= 0 or ny <= 0 or nz <= 0: error_all("Invalid grid size")

  if inpx == 0 and inpy == 0 and inpz == 0: pass
  elif inpx <= 0 or inpy <= 0 or inpz <= 0: error_all("Invalid proc grid")
  elif inpx*inpy != nprocs: error_all("Specified proc grid does not match nprocs")

  if outpx == 0 and outpy == 0 and outpz == 0: pass
  elif outpx <= 0 or outpy <= 0 or outpy <= 0: error_all("Invalid proc grid")
  elif outpx*outpy != nprocs:
    error_all("Specified proc grid does not match nprocs")

  if nloop <= 0: error_all("Invalid Nloop")
  if nloop == 0 and tuneflag == 0: error_all("Invalid Nloop")
  if iflag == RANDOM and seed <= 0: error_all("Invalid initialize setting")
  if mode < 0 or mode > 3: error_all("Invalid FFT mode")
  if mode > 1 and vflag: error_all("Cannot validate forward only FFT")

  if tuneflag and tuneper <= 0: error_all("Invalid tune nper")
  if tuneflag and tunemax < 0.0: error_all("Invalid tune tmax")
  if tuneflag and (tuneextra < 0 or tuneextra > 1):
    error_all("Invalid tune extra")
  if tuneflag and rflag: error_all("Cannot tune with remap only")

# ----------------------------------------------------------------------
# partition processors across grid dimensions
# flag = IN for input partitions, or OUT for output partitions
# if user set Px,Py,Pz -> just return
# for IN:
#   assign nprocs as bricks to 3d grid to minimize surface area per proc
#   derived from SPPARKS Domain::procs2domain_3d()
# for OUT:
#   assign nprocs as rectangles to xy grid to minimize surface area per proc
#   derived from SPPARKS Domain::procs2domain_2d()
# ----------------------------------------------------------------------

def proc_setup(flag):
  global inpx,inpy,inpz,outpx,outpy,outpz

  if flag == IN:
    if inpx or inpy or inpz: return
    inpx,inpy,inpz = proc3d()
    
  if flag == OUT:
    if outpx or outpy or outpz: return
    if mode == 0 or mode == 2: outpx,outpy,outpz = proc3d()
    if mode == 1 or mode == 3: outpx,outpy,outpz = proc2d()

def proc3d():
  xprd = nx
  yprd = ny
  zprd = nz
  bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd)
  
  # loop thru all possible factorizations of nprocs
  # surf = surface area of a proc sub-domain
  
  ipx = 1
  while ipx <= nprocs:
    if nprocs % ipx == 0:
      nremain = nprocs/ipx
      ipy = 1
      while ipy <= nremain:
        if nremain % ipy == 0:
          ipz = nremain/ipy
          boxx = xprd/ipx
          boxy = yprd/ipy
          boxz = zprd/ipz
          surf = boxx*boxy + boxy*boxz + boxz*boxx
          if surf < bestsurf:
            bestsurf = surf
            px = ipx
            py = ipy
            pz = ipz
        ipy += 1
    ipx += 1

  if px*py*pz != nprocs: error_all("Computed proc grid does not match nprocs")
  return px,py,pz

def proc2d():
  xprd = nx
  yprd = ny
  bestsurf = 2.0 * (xprd+yprd)
  
  # loop thru all possible factorizations of nprocs
  # surf = surface area of a proc sub-domain
  
  ipx = 1
  while ipx <= nprocs:
    if nprocs % ipx == 0:
      ipy = nprocs/ipx
      boxx = xprd/ipx
      boxy = yprd/ipy
      surf = boxx + boxy
      if surf < bestsurf:
        bestsurf = surf
        px = ipx
        py = ipy
    ipx += 1

  pz = 1
  if px*py*pz != nprocs: error_all("Computed proc grid does not match nprocs")
  return px,py,pz

# ----------------------------------------------------------------------
# partition FFT grid
# once for input grid, once for output grid
# use Px,Py,Pz for in/out
# ----------------------------------------------------------------------

def grid_setup():
  global inxlo,inxhi,inylo,inyhi,inzlo,inzhi
  global outxlo,outxhi,outylo,outyhi,outzlo,outzhi
  global nfft_in,nfft_out
   
  # ipx,ipy,ipz = my position in input 2d grid of procs

  ipx = me % inpx
  ipy = (me/inpx) % inpy
  ipz = me / (inpx*inpy)

  # nlo,nhi = lower/upper limits of the 2d brick I own

  inxlo = int(1.0 * ipx * nx / inpx)
  inxhi = int(1.0 * (ipx+1) * nx / inpx) - 1

  inylo = int(1.0 * ipy * ny / inpy)
  inyhi = int(1.0 * (ipy+1) * ny / inpy) - 1

  inzlo = int(1.0 * ipz * nz / inpz)
  inzhi = int(1.0 * (ipz+1) * nz / inpz) - 1

  nfft_in = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1)

  # ipx,ipy = my position in output 2d grid of procs

  ipx = me % outpx
  ipy = (me/outpx) % outpy
  ipz = me / (outpx*outpy)

  # nlo,nhi = lower/upper limits of the 2d brick I own

  outxlo = int(1.0 * ipx * nx / outpx)
  outxhi = int(1.0 * (ipx+1) * nx / outpx) - 1

  outylo = int(1.0 * ipy * ny / outpy)
  outyhi = int(1.0 * (ipy+1) * ny / outpy) - 1

  outzlo = int(1.0 * ipz * nz / outpz)
  outzhi = int(1.0 * (ipz+1) * nz / outpz) - 1

  nfft_out = (outxhi-outxlo+1) * (outyhi-outylo+1) * (outzhi-outzlo+1)

# ----------------------------------------------------------------------
# create FFT plan
# ----------------------------------------------------------------------

def plan():
  global fft,fftsize,sendsize,recvsize,timesetup,timetune,nloop
  
  fft = FFT3dMPI(world,precision)
  fft.set("remaponly",rflag)

  fft.set("collective",cflag)
  fft.set("exchange",eflag)
  fft.set("pack",pflag)

  if mode == 0 or mode == 2: permute = 0
  else: permute = 2

  # will use fftsize to allocate work buffer
  # ignore sendsize, recvsize b/c let FFT allocate remap buffers internally
  # set timesetup and timetune
  # reset nloop if tuning and user nloop = 0

  world.Barrier()
  time1 = MPI.Wtime()

  if not tuneflag:
    fftsize,sendsize,recvsize = \
      fft.setup(nx,ny,nz,inxlo,inxhi,inylo,inyhi,inzlo,inzhi,
                outxlo,outxhi,outylo,outyhi,outzlo,outzhi,permute)
  else:
    flag = 0
    if mode >= 2: flag = 1
    fftsize,sendsize,recvsize = fft.tune(nx,ny,nz,
              inxlo,inxhi,inylo,inyhi,inzlo,inzhi,
              outxlo,outxhi,outylo,outyhi,outzlo,outzhi,
              permute,flag,tuneper,tunemax,tuneextra)
    if nloop == 0: nloop = fft.get_int("npertrial")

  world.Barrier()
  time2 = MPI.Wtime()

  if not tuneflag:
    timesetup = time2 - time1
    timetune = 0.0
  else:
    timesetup = fft.get_double("setuptime")
    timetune = time2 - time1

# ----------------------------------------------------------------------
# allocate memory for FFT grid
# ----------------------------------------------------------------------

def allocate():
  global work
  if precision == 1: work = np.zeros(2*fftsize,np.float32)
  else: work = np.zeros(2*fftsize,np.float)

# ----------------------------------------------------------------------
# initialize FFT grid
# ----------------------------------------------------------------------

def initialize():
  if iflag == ZERO:
    for m in xrange(2*nfft_in): work[m] = 0.0

  elif iflag == STEP:
    nxlocal = inxhi - inxlo + 1
    nylocal = inyhi - inylo + 1

    for m in xrange(nfft_in):
      ilocal = m % nxlocal
      jlocal = (m/nxlocal) % nylocal
      klocal = m / (nxlocal*nylocal)
      iglobal = inxlo + ilocal
      jglobal = inylo + jlocal
      kglobal = inzlo + klocal
      if iglobal < nx/2 and jglobal < ny/2 and kglobal < nz/2: work[2*m] = 1.0
      else: work[2*m] = 0.0
      work[2*m+1] = 0.0

  elif iflag == INDEX:
    nxlocal = inxhi - inxlo + 1
    nylocal = inyhi - inylo + 1

    for m in xrange(nfft_in):
      ilocal = m % nxlocal
      jlocal = (m/nxlocal) % nylocal
      klocal = m / (nxlocal*nylocal)
      iglobal = inxlo + ilocal
      jglobal = inylo + jlocal
      kglobal = inzlo + klocal
      work[2*m] = kglobal + jglobal + iglobal + 1
      work[2*m+1] = 0.0

  elif iflag == RANDOM:
    for m in xrange(2*nfft_in):
      work[m] = random()

# ----------------------------------------------------------------------
# output FFT grid values
# flag = 0 for initial partition
# flag = 1 for final partition
# ----------------------------------------------------------------------

def output(flag,txt):

  if me == 0: print txt

  for iproc in range(nprocs):
    if me != iproc: continue
    if me >= 1: tmp = world.Recv(source=me-1,tag=0)

    if flag == 0:
      nxlocal = inxhi - inxlo + 1
      for m in xrange(nfft_in):
        ilocal = m % nxlocal
        jlocal = m / nxlocal
        iglobal = inxlo + ilocal
        jglobal = inylo + jlocal
        print "Value (%d,%d) on proc %d = (%g,%g)" % \
          (iglobal,jglobal,me,work[2*m],work[2*m+1])
    else:
      nxlocal = outxhi - outxlo + 1
      for m in xrange(nfft_in):
        ilocal = m % nxlocal
        jlocal = m / nxlocal
        iglobal = outxlo + ilocal
        jglobal = outylo + jlocal
        print "Value (%d,%d) on proc %d = (%g,%g)" % \
          (iglobal,jglobal,me,work[2*m],work[2*m+1])

    if me < nprocs-1: world.send(tmp,dest=me+1,tag=0)

# ----------------------------------------------------------------------
# validation check for correct result
# ----------------------------------------------------------------------

def validate():
  global seed,epsmax

  epsilon = 0.0
  
  if iflag == ZERO:
    for m in xrange(2*nfft_in):
      delta = math.fabs(work[m])
      if delta > epsilon: epsilon = delta

  elif iflag == STEP:
    nxlocal = inxhi - inxlo + 1
    nylocal = inyhi - inylo + 1
    
    for m in xrange(nfft_in):
      ilocal = m % nxlocal
      jlocal = (m/nxlocal) % nylocal
      klocal = m / (nxlocal*nylocal)
      iglobal = inxlo + ilocal
      jglobal = inylo + jlocal
      kglobal = inzlo + klocal
      if iglobal < nx/2 and jglobal < ny/2 and kglobal < nz/2: value = 1.0
      else: value = 0.0
      delta = math.fabs(work[2*m]-value)
      if delta > epsilon: epsilon = delta
      delta = math.fabs(work[2*m+1])
      if delta > epsilon: epsilon = delta

  elif iflag == INDEX:
    nxlocal = inxhi - inxlo + 1
    nylocal = inyhi - inylo + 1
    
    for m in xrange(nfft_in):
      ilocal = m % nxlocal
      jlocal = (m/nxlocal) % nylocal
      klocal = m / (nxlocal*nylocal)
      iglobal = inxlo + ilocal
      jglobal = inylo + jlocal
      kglobal = inzlo + klocal
      value = kglobal + jglobal + iglobal + 1
      delta = math.fabs(work[2*m]-value)
      if delta > epsilon: epsilon = delta
      delta = math.fabs(work[2*m+1])
      if delta > epsilon: epsilon = delta

  elif iflag == RANDOM:
    seed = seedinit
    for m in xrange(2*nfft_in):
      newvalue = random()
      delta = math.fabs(work[m]-newvalue)
      if delta > epsilon: epsilon = delta
      
  epsmax = world.allreduce(epsilon,op=MPI.MAX)

# ----------------------------------------------------------------------
# output timing data
# ----------------------------------------------------------------------

def timing():

  # perform only 1d FFTs

  if tflag:
    for i in xrange(2*nfft_in): work[i] = 0.0

    world.Barrier()
    time1 = MPI.Wtime()

    if mode < 2:
      for i in xrange(nloop):
        fft.only_1d_ffts(work,1)
        fft.only_1d_ffts(work,-1)
    else:
      for i in xrange(nloop):
        fft.only_1d_ffts(work,1)

    world.Barrier()
    time2 = MPI.Wtime()
    time1d = time2 - time1

  # perform all remaps

  if tflag:
    for i in xrange(2*nfft_in): work[i] = 0.0

    world.Barrier()
    time1 = MPI.Wtime()

    if mode < 2:
      for i in xrange(nloop):
        fft.only_remaps(work,work,1)
        fft.only_remaps(work,work,-1)
    else:
      for i in xrange(nloop):
        fft.only_remaps(work,work,1)

    world.Barrier()
    time2 = MPI.Wtime()
    time_remap = time2 - time1

  # perform only single remaps

  if tflag:
    for i in xrange(2*nfft_in): work[i] = 0.0

    world.Barrier()
    time1 = MPI.Wtime()

    if mode < 2:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,1)
        fft.only_one_remap(work,work,-1,1)
    else:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,1)

    world.Barrier()
    time2 = MPI.Wtime()
    time_remap1 = time2 - time1

    if mode < 2:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,2)
        fft.only_one_remap(work,work,-1,2)
    else:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,2)

    world.Barrier()
    time3 = MPI.Wtime()
    time_remap2 = time3 - time2

    if mode < 2:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,3)
        fft.only_one_remap(work,work,-1,3)
    else:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,3)

    world.Barrier()
    time4 = MPI.Wtime()
    time_remap3 = time4 - time3

    if mode < 2:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,4)
        fft.only_one_remap(work,work,-1,4)
    else:
      for i in xrange(nloop):
        fft.only_one_remap(work,work,1,4)

    world.Barrier()
    time5 = MPI.Wtime()
    time_remap4 = time5 - time4

  # stats output
  # nfft = 2x larger for modes 0,1

  if mode < 2: nfft = 2*nloop
  else: nfft = nloop

  onetime = timefft/nfft
  nsize = nx * ny * nz
  log2n = math.log(nsize)/math.log(2.0)
  floprate = 5.0 * nsize * log2n / onetime / (1024*1024*1024)
  gridbytes = 4*precision * 2*fftsize
  
  if me == 0:
    print "3d FFTs with %s library, precision = %s" % \
      (fft.get_string("fft1d"),fft.get_string("precision"))
    print "Grid size: %d %d %d" % (nx,ny,nz)
    print "  initial proc grid: %d %d %d" % (inpx,inpy,inpz)
    print "  x pencil proc grid: %d %d %d" % \
       (fft.get_int("npfast1"),fft.get_int("npfast2"),fft.get_int("npfast3"))
    print "  y pencil proc grid: %d %d %d" % \
       (fft.get_int("npmid1"),fft.get_int("npmid2"),fft.get_int("npmid3"))
    print "  z pencil proc grid: %d %d %d" % \
       (fft.get_int("npslow1"),fft.get_int("npslow2"),fft.get_int("npslow3"))
    print "  2d brick proc grid: %d %d %d" % \
       (fft.get_int("npbrick1"),fft.get_int("npbrick2"),fft.get_int("npbrick3"))
    print "  final proc grid: %d %d %d" % (outpx,outpy,outpz)
    
    if tuneflag:
      ntrial = fft.get_int("ntrial")
      npertrial = fft.get_int("npertrial")
      print "Tuning trials & iterations: %d %d" % (ntrial,npertrial)
      for i in range(ntrial):
        print "  coll exch pack 3dFFT 1dFFT remap r1 r2 r3 r4: " + \
          "%d %d %d %g %g %g %g %g %g %g" % \
          (fft.get_int_vector("cflags")[i],fft.get_int_vector("eflags")[i],
           fft.get_int_vector("pflags")[i],fft.get_double_vector("tfft")[i],
           fft.get_double_vector("t1d")[i],fft.get_double_vector("tremap")[i],
           fft.get_double_vector("tremap1")[i],
           fft.get_double_vector("tremap2")[i],
           fft.get_double_vector("tremap3")[i],
           fft.get_double_vector("tremap4")[i])

    if mode == 0:
      print "%d forward and %d back FFTs on %d procs" % (nloop,nloop,nprocs)
    elif mode == 1:
      print "%d forward and %d back convolution FFTs on %d procs" % \
        (nloop,nloop,nprocs)
    elif mode == 2:
      print "%d forward FFTs on %d procs" % (nloop,nprocs)
    elif mode == 3:
      print "%d forward convolution FFTs on %d procs" % (nloop,nprocs)

    print "Collective, exchange, pack methods: %d %d %d" % \
      (fft.get_int("collective"),fft.get_int("exchange"),fft.get_int("pack"))
    print "Memory usage (per-proc) for FFT grid = %g MBytes" % \
      (float(gridbytes) / 1024/1024)
    print "Memory usage (per-proc) by fftMPI = %g MBytes" % \
      (float(fft.get_int64("memusage")) / 1024/1024)

    if vflag: print "Max error = %g" % epsmax
    if not tuneflag: print "Initialize grid = %g secs" % (timeinit-timesetup)
    else: print "Initialize grid = %g secs" % (timeinit-timetune)
    print "FFT setup = %g secs" % timesetup
    print "FFT tune = %g secs" % timetune
    print "Time for 3d FFTs = %g secs" % timefft
    print "  time/fft2d = %g secs" % onetime
    print "  flop rate for 3d FFTs = %g Gflops" % floprate
    if tflag:
      print "Time for 1d FFTs only = %g secs" % time1d
      print "  time/fft1d = %g secs" % (time1d/nfft)
      print "  fraction of time in 1d FFTs = %g" % (time1d/timefft)
    if tflag:
      print "Time for remaps only = %g secs" % time_remap
      print "  fraction of time in remaps = %g" % (time_remap/timefft)
      print "Time for remap #1 = %g secs" % time_remap1
      print "  fraction of time in remap #1 = %g" % (time_remap1/timefft)
      print "Time for remap #2 = %g secs" % time_remap2
      print "  fraction of time in remap #2 = %g" % (time_remap2/timefft)
      print "Time for remap #3 = %g secs" % time_remap3
      print "  fraction of time in remap #3 = %g" % (time_remap3/timefft)
      print "Time for remap #4 = %g secs" % time_remap4
      print "  fraction of time in remap #4 = %g" % (time_remap4/timefft)

# ----------------------------------------------------------------------
# deallocate memory for FFT grid
# no need to do this since Python does garbage collection
# ----------------------------------------------------------------------

def deallocate():
  pass

# ----------------------------------------------------------------------
# main program
# ----------------------------------------------------------------------
  
# MPI setup

world = MPI.COMM_WORLD
me = world.rank
nprocs = world.size

# parse command-line args

options()

# partition FFT grid across procs, for both input and output
# create FFT plan, tune if requested
# allocate grid
# initialize FFT grid
# grid output

world.Barrier()
time1 = MPI.Wtime()

proc_setup(IN)
proc_setup(OUT)
grid_setup()
plan()
allocate()
initialize()

world.Barrier()
time2 = MPI.Wtime()
timeinit = time2 - time1

if oflag: output(0,"Initial grid")

# perform FFTs

world.Barrier()
time1 = MPI.Wtime()

if mode < 2:
  for i in range(nloop):
    fft.compute(work,work,1)
    # if oflag: output(1,"Middle grid")
    fft.compute(work,work,-1)
else:
  for i in range(nloop):
    fft3d_compute(fft,work,work,1)

world.Barrier()
time2 = MPI.Wtime()
timefft = time2 - time1

# validation check on result
# grid output
# timing results
# deallocate grid and plan

if vflag: validate()
if oflag:
  if mode < 2: output(0,"Final grid")
  else: output(1,"Final grid")
timing()
deallocate()
del fft
