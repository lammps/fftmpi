#!/usr/bin/env python

# Compute a forward/inverse complex FFT using fftMPI
#   change FFT size by editing 3 "FFT size" lines
#   run on any number of procs

# Run syntax:
# % python simple,py                 # run in serial
# % mpirun -np 4 pythone simple.py   # run in parallel

import math
import numpy as np
from mpi4py import MPI
from fftmpi import FFT3dMPI

# FFT size

NFAST = 128
NMID = 128
NSLOW = 128

# precision-dependent settings
# NOTE: how to set this at run-time

precision = 2

# ----------------------------------------------------------------------
# main program
# ----------------------------------------------------------------------

# setup MPI

world = MPI.COMM_WORLD
me = world.rank
nprocs = world.size

# instantiate FFT

fft = FFT3dMPI(world,precision)

# simple algorithm to factor Nprocs into roughly cube roots

npfast = int(math.pow(nprocs,1.0/3.0))
while npfast < nprocs:
  if nprocs % npfast == 0: break
  npfast += 1
npmidslow = nprocs / npfast
npmid = int(math.sqrt(npmidslow))
while npmid < npmidslow:
  if npmidslow % npmid == 0: break
  npmid += 1
npslow = nprocs / npfast / npmid

# partition grid into Npfast x Npmid x Npslow bricks

nfast = NFAST
nmid = NMID
nslow = NSLOW

ipfast = me % npfast
ipmid = (me/npfast) % npmid
ipslow = me / (npfast*npmid)

ilo = int(1.0*ipfast*nfast/npfast)
ihi = int(1.0*(ipfast+1)*nfast/npfast) - 1
jlo = int(1.0*ipmid*nmid/npmid)
jhi = int(1.0*(ipmid+1)*nmid/npmid) - 1
klo = int(1.0*ipslow*nslow/npslow)
khi = int(1.0*(ipslow+1)*nslow/npslow) - 1

# setup FFT, could replace with tune()

fftsize,sendsize,recvsize = \
  fft.setup(nfast,nmid,nslow,
            ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi,0)

# tune FFT, could replace with setup()

#fftsize,sendsize,recvsize = \
#  fft.tune(nfast,nmid,nslow,
#           ilo,ihi,jlo,jhi,klo,khi,ilo,ihi,jlo,jhi,klo,khi,0,0,5,10.0,0)

# initialize each proc's local grid
# global initialization is specific to proc count

if precision == 1: work = np.zeros(2*fftsize,np.float32)
else: work = np.zeros(2*fftsize,np.float)

n = 0
for k in xrange(klo,khi+1):
  for j in xrange(jlo,jhi+1):
    for i in xrange(ilo,ihi+1):
      work[n] = n
      n += 1
      work[n] = n
      n += 1

# perform 2 FFTs

timestart = MPI.Wtime()
fft.compute(work,work,1)
fft.compute(work,work,-1)
timestop = MPI.Wtime()

if me == 0:
  print "Two %dx%dx%d FFTs on %d procs as %dx%dx%d grid" % \
    (nfast,nmid,nslow,nprocs,npfast,npmid,npslow)
  print "CPU time = %g secs" % (timestop-timestart)

# find largest difference between initial/final values
# should be near zero

n = 0
mydiff = 0.0
for k in xrange(klo,khi+1):
  for j in xrange(jlo,jhi+1):
    for i in xrange(ilo,ihi+1):
      if abs(work[n]-n) > mydiff: mydiff = abs(work[n]-n)
      n += 1
      if abs(work[n]-n) > mydiff: mydiff = abs(work[n]-n)
      n += 1

world.allreduce(mydiff,op=MPI.MAX)
if me == 0: print "Max difference in initial/final values =",mydiff

# clean up

del fft
