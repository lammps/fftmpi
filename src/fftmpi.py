# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

# Python wrapper on fftMPI library via ctypes

import sys,traceback
from ctypes import *

# Numpy and mpi4py packages must exist
  
try:
  import numpy as np
  numpyflag = 1
except:
  print "fftMPI: Must have Numpy installed in python to use fftMPI"
  sys.exit()
  
try:
  from mpi4py import MPI
  mpi4pyflag = 1
except:
  print "fftMPI: Cannot pass MPI communicator w/out mpi4py package"
  sys.exit()

# ------------------------------------------------------------------
# instantiate fft3dMPI thru its C-interface
# ------------------------------------------------------------------

class FFT3dMPI:
  
  def __init__(self,comm,precision):

    self.precision = precision
    if precision != 1 and precision != 2:
      # NOTE: how to handle errors, e.g. in parallel
      print "fftMPI: precision must be 1 or 2"
      sys.exit()
    
    # load fft3dmpi.so
    
    try:
      self.lib = CDLL("libfft3dmpi.so",RTLD_GLOBAL)
    except:
      etype,value,tb = sys.exc_info()
      traceback.print_exception(etype,value,tb)
      raise OSError,"Could not load fft3dMPI dynamic library"

    # define ctypes API for each library method

    if mpi4pyflag and MPI._sizeof(comm) == sizeof(c_int): MPI_Comm = c_int
    else: MPI_Comm = c_void_p
    
    self.lib.fft3d_create.argtypes = [MPI_Comm,c_int,POINTER(c_void_p)]
    self.lib.fft3d_create.restype = None

    self.lib.fft3d_destroy.argtypes = [c_void_p]
    self.lib.fft3d_destroy.restype = None

    self.lib.fft3d_set.argtypes = [c_void_p,c_char_p,c_int]
    self.lib.fft3d_set.restype = None
    
    self.lib.fft3d_get_int.argtypes = [c_void_p,c_char_p]
    self.lib.fft3d_get_int.restype = c_int
    self.lib.fft3d_get_int64.argtypes = [c_void_p,c_char_p]
    self.lib.fft3d_get_int64.restype = c_longlong
    self.lib.fft3d_get_double.argtypes = [c_void_p,c_char_p]
    self.lib.fft3d_get_double.restype = c_double
    self.lib.fft3d_get_string.argtypes = [c_void_p,c_char_p,POINTER(c_int)]
    self.lib.fft3d_get_string.restype = c_char_p
    self.lib.fft3d_get_int_vector.argtypes = [c_void_p,c_char_p,POINTER(c_int)]
    self.lib.fft3d_get_int_vector.restype = POINTER(c_int)
    self.lib.fft3d_get_double_vector.argtypes = [c_void_p,c_char_p,POINTER(c_int)]
    self.lib.fft3d_get_double_vector.restype = POINTER(c_double)

    self.lib.fft3d_setup.argtypes = \
      [c_void_p,c_int,c_int,c_int,
       c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,POINTER(c_int),POINTER(c_int),
       POINTER(c_int)]
    self.lib.fft3d_setup.restype = None

    if precision == 1:
      self.lib.fft3d_setup_memory.argtypes = \
        [c_void_p,POINTER(c_float),
         POINTER(c_float)]
    elif precision == 2:
      self.lib.fft3d_setup_memory.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double)]
    self.lib.fft3d_setup_memory.restype = None

    if precision == 1:
      self.lib.fft3d_compute.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),c_int]
    elif precision == 2:
      self.lib.fft3d_compute.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),c_int]
    self.lib.fft3d_compute.restype = None

    if precision == 1:
      self.lib.fft3d_only_1d_ffts.argtypes = \
        [c_void_p,POINTER(c_float),c_int]
    elif precision == 2:
      self.lib.fft3d_only_1d_ffts.argtypes = \
        [c_void_p,POINTER(c_double),c_int]
    self.lib.fft3d_only_1d_ffts.restype = None

    if precision == 1:
      self.lib.fft3d_only_remaps.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),c_int]
    elif precision == 2:
      self.lib.fft3d_only_remaps.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),c_int]
    self.lib.fft3d_only_remaps.restype = None

    if precision == 1:
      self.lib.fft3d_only_one_remap.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),c_int,c_int]
    elif precision == 2:
      self.lib.fft3d_only_one_remap.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),c_int,c_int]
    self.lib.fft3d_only_remaps.restype = None

    self.lib.fft3d_tune.argtypes = \
      [c_void_p,c_int,c_int,c_int,
       c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,POINTER(c_int),POINTER(c_int),POINTER(c_int),
       c_int,c_int,c_double,c_int]
    self.lib.fft3d_tune.restype = None

    # create an instance of fftMPI

    self.fft = c_void_p()
    comm_ptr = MPI._addressof(comm)
    comm_value = MPI_Comm.from_address(comm_ptr)
    self.lib.fft3d_create(comm_value,precision,byref(self.fft))

  # destroy instance of fftMPI
  
  def __del__(self):
    self.lib.fft3d_destroy(self.fft)

  def destroy(self):
    self.lib.fft3d_destroy(self.fft)
    self.lib = None

  # setup

  def set(self,keyword,value):
    self.lib.fft3d_set(self.fft,keyword,value)

  def get_int(self,keyword):
    return self.lib.fft3d_get_int(self.fft,keyword)

  def get_int64(self,keyword):
    return self.lib.fft3d_get_int64(self.fft,keyword)

  def get_double(self,keyword):
    return self.lib.fft3d_get_double(self.fft,keyword)

  def get_string(self,keyword):
    nlen = c_int()
    cptr = self.lib.fft3d_get_string(self.fft,keyword,byref(nlen))
    if not bool(cptr): return None   # NULL ptr is False in ctypes
    return cptr                      # becomes Python string with correct length
    
  def get_int_vector(self,keyword):
    nlen = c_int()
    cptr = self.lib.fft3d_get_int_vector(self.fft,keyword,byref(nlen))
    if not bool(cptr): return None   # NULL ptr is False in ctypes
    ivec = cptr[:nlen.value]         # ivec now has correct Python length
    return ivec
    
  def get_double_vector(self,keyword):
    nlen = c_int()
    cptr = self.lib.fft3d_get_double_vector(self.fft,keyword,byref(nlen))
    if not bool(cptr): return None   # NULL ptr is False in ctypes
    dvec = cptr[:nlen.value]         # dvec now has correct Python length
    return dvec

  def setup(self,nfast,nmid,nslow,
            in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
            out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,permute):
    fftsize = c_int()
    sendsize = c_int()
    recvsize = c_int()
    self.lib.fft3d_setup(self.fft,nfast,nmid,nslow,
                         in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                         out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                         permute,byref(fftsize),byref(sendsize),byref(recvsize))
    fftsize = fftsize.value
    sendsize = sendsize.value
    recvsize = recvsize.value
    return fftsize,sendsize,recvsize

  def setup_memory(self,sendbuf,recvbuf):
    if "numpy" not in str(type(sendbuf)) or "numpy" not in str(type(recvbuf)):
      print "Must use Numpy arrays for setup_memory() data"
      sys.exit()
      
    if self.precision == 1:
      csendbuf = csendbuf.ctypes.data_as(POINTER(c_float))
      crecvbuf = crecvbuf.ctypes.data_as(POINTER(c_float))
    else:
      csendbuf = csendbuf.ctypes.data_as(POINTER(c_double))
      crecvbuf = crecvbuf.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft3d_setup_memory(self.fft,csendbuf,crecvbuf)
      
  # compute

  def compute(self,indata,outdata,flag):
    if "numpy" not in str(type(indata)) or "numpy" not in str(type(outdata)):
      print "Must use Numpy arrays for compute() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft3d_compute(self.fft,cin,cout,flag)

  def only_1d_ffts(self,indata,flag):
    if "numpy" not in str(type(indata)):
      print "Must use Numpy arrays for only_1d_ffts() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft3d_only_1d_ffts(self.fft,cin,flag)

  def only_remaps(self,indata,outdata,flag):
    if "numpy" not in str(type(indata)) or "numpy" not in str(type(outdata)):
      print "Must use Numpy arrays for only_remaps() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))

    self.lib.fft3d_only_remaps(self.fft,cin,cout,flag)

  def only_one_remap(self,indata,outdata,flag,which):
    if "numpy" not in str(type(indata)) or "numpy" not in str(type(outdata)):
      print "Must use Numpy arrays for only_one_remap() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft3d_only_one_remap(self.fft,cin,cout,flag,which)

  # tune

  def tune(self,nfast,nmid,nslow,
           in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
           out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
           permute,flag,niter,tmax,tflag):
    fftsize = c_int()
    sendsize = c_int()
    recvsize = c_int()
    self.lib.fft3d_tune(self.fft,nfast,nmid,nslow,
                        in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                        out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                        permute,byref(fftsize),byref(sendsize),byref(recvsize),
                        flag,niter,tmax,tflag)
    fftsize = fftsize.value
    sendsize = sendsize.value
    recvsize = recvsize.value
    return fftsize,sendsize,recvsize

# ------------------------------------------------------------------
# instantiate fft2dMPI thru its C-interface
# ------------------------------------------------------------------

class FFT2dMPI:

  def __init__(self,comm,precision):

    self.precision = precision
    if precision != 1 and precision != 2:
      # NOTE: how to handle errors, e.g. in parallel
      print "fftMPI: precision must be 1 or 2"
      sys.exit()
    
    # load fft2dmpi.so
    
    try:
      self.lib = CDLL("libfft2dmpi.so",RTLD_GLOBAL)
    except:
      etype,value,tb = sys.exc_info()
      traceback.print_exception(etype,value,tb)
      raise OSError,"Could not load fft2dMPI dynamic library"

    # define ctypes API for each library method

    if mpi4pyflag and MPI._sizeof(comm) == sizeof(c_int): MPI_Comm = c_int
    else: MPI_Comm = c_void_p

    self.lib.fft2d_create.argtypes = [MPI_Comm,c_int,POINTER(c_void_p)]
    self.lib.fft2d_create.restype = None

    self.lib.fft2d_destroy.argtypes = [c_void_p]
    self.lib.fft2d_destroy.restype = None

    self.lib.fft2d_set.argtypes = [c_void_p,c_char_p,c_int]
    self.lib.fft2d_set.restype = None

    self.lib.fft2d_get_int.argtypes = [c_void_p,c_char_p]
    self.lib.fft2d_get_int.restype = c_int
    self.lib.fft2d_get_int64.argtypes = [c_void_p,c_char_p]
    self.lib.fft2d_get_int64.restype = c_longlong
    self.lib.fft2d_get_double.argtypes = [c_void_p,c_char_p]
    self.lib.fft2d_get_double.restype = c_double
    self.lib.fft2d_get_string.argtypes = [c_void_p,c_char_p,POINTER(c_int)]
    self.lib.fft2d_get_string.restype = c_char_p
    self.lib.fft2d_get_int_vector.argtypes = [c_void_p,c_char_p,POINTER(c_int)]
    self.lib.fft2d_get_int_vector.restype = POINTER(c_int)
    self.lib.fft2d_get_double_vector.argtypes = [c_void_p,c_char_p,POINTER(c_int)]
    self.lib.fft2d_get_double_vector.restype = POINTER(c_double)
    
    self.lib.fft2d_setup.argtypes = \
      [c_void_p,c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,POINTER(c_int),POINTER(c_int),POINTER(c_int)]
    self.lib.fft2d_setup.restype = None

    if precision == 1:
      self.lib.fft2d_setup_memory.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float)]
    elif precision == 2:
      self.lib.fft2d_setup_memory.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double)]
    self.lib.fft2d_setup_memory.restype = None

    if precision == 1:
      self.lib.fft2d_compute.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),c_int]
    elif precision == 2:
      self.lib.fft2d_compute.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),c_int]
    self.lib.fft2d_compute.restype = None

    if precision == 1:
      self.lib.fft2d_only_1d_ffts.argtypes = \
        [c_void_p,POINTER(c_float),c_int]
    elif precision == 2:
      self.lib.fft2d_only_1d_ffts.argtypes = \
        [c_void_p,POINTER(c_double),c_int]
    self.lib.fft2d_only_1d_ffts.restype = None

    if precision == 1:
      self.lib.fft2d_only_remaps.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),c_int]
    elif precision == 2:
      self.lib.fft2d_only_remaps.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),c_int]
    self.lib.fft2d_only_remaps.restype = None

    if precision == 1:
      self.lib.fft2d_only_one_remap.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),c_int,c_int]
    elif precision == 2:
      self.lib.fft2d_only_one_remap.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),c_int,c_int]
    self.lib.fft2d_only_remaps.restype = None

    self.lib.fft2d_tune.argtypes = \
    [c_void_p,c_int,c_int,
       c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,POINTER(c_int),POINTER(c_int),POINTER(c_int),
       c_int,c_int,c_double,c_int]
    self.lib.fft2d_tune.restype = None

    # create an instance of fftMPI

    self.fft = c_void_p()
    comm_ptr = MPI._addressof(comm)
    comm_value = MPI_Comm.from_address(comm_ptr)
    self.lib.fft2d_create(comm_value,precision,byref(self.fft))

  # destroy instance of fftMPI
  
  def __del__(self):
    self.lib.fft2d_destroy(self.fft)

  def destroy(self):
    self.lib.fft2d_destroy(self.fft)
    self.lib = None

  # setup

  def set(self,keyword,value):
    self.lib.fft2d_set(self.fft,keyword,value)

  def get_int(self,keyword):
    return self.lib.fft2d_get_int(self.fft,keyword)

  def get_int64(self,keyword):
    return self.lib.fft2d_get_int64(self.fft,keyword)

  def get_double(self,keyword):
    return self.lib.fft2d_get_double(self.fft,keyword)

  def get_string(self,keyword):
    nlen = c_int()
    cptr = self.lib.fft2d_get_string(self.fft,keyword,byref(nlen))
    if not bool(cptr): return None   # NULL ptr is False in ctypes
    return cptr                      # becomes Python string with correct length
    
  def get_int_vector(self,keyword):
    nlen = c_int()
    cptr = self.lib.fft2d_get_int_vector(self.fft,keyword,byref(nlen))
    if not bool(cptr): return None   # NULL ptr is False in ctypes
    ivec = cptr[:nlen.value]         # ivec now has correct Python length
    return ivec
    
  def get_double_vector(self,keyword):
    nlen = c_int()
    cptr = self.lib.fft2d_get_double_vector(self.fft,keyword,byref(nlen))
    if not bool(cptr): return None   # NULL ptr is False in ctypes
    dvec = cptr[:nlen.value]         # dvec now has correct Python length
    return dvec

  def setup(self,nfast,nslow,in_ilo,in_ihi,in_jlo,in_jhi,
            out_ilo,out_ihi,out_jlo,out_jhi,permute):
    fftsize = c_int()
    sendsize = c_int()
    recvsize = c_int()
    self.lib.fft2d_setup(self.fft,nfast,nslow,in_ilo,in_ihi,in_jlo,in_jhi,
                         out_ilo,out_ihi,out_jlo,out_jhi,
                         permute,byref(fftsize),byref(sendsize),byref(recvsize))
    fftsize = fftsize.value
    sendsize = sendsize.value
    recvsize = recvsize.value
    return fftsize,sendsize,recvsize

  def setup_memory(self,sendbuf,recvbuf):
    if "numpy" not in str(type(sendbuf)) or "numpy" not in str(type(recvbuf)):
      print "Must use Numpy arrays for setup_memory() data"
      sys.exit()
      
    if self.precision == 1:
      csendbuf = csendbuf.ctypes.data_as(POINTER(c_float))
      crecvbuf = crecvbuf.ctypes.data_as(POINTER(c_float))
    else:
      csendbuf = csendbuf.ctypes.data_as(POINTER(c_double))
      crecvbuf = crecvbuf.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft2d_setup_memory(self.fft,csendbuf,crecvbuf)
      
  # compute

  def compute(self,indata,outdata,flag):
    if "numpy" not in str(type(indata)) or "numpy" not in str(type(outdata)):
      print "Must use Numpy arrays for compute() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft2d_compute(self.fft,cin,cout,flag)

  def only_1d_ffts(self,indata,flag):
    if "numpy" not in str(type(indata)):
      print "Must use Numpy arrays for only_1d_ffts() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft2d_only_1d_ffts(self.fft,cin,flag)

  def only_remaps(self,indata,outdata,flag):
    if "numpy" not in str(type(indata)) or "numpy" not in str(type(outdata)):
      print "Must use Numpy arrays for only_remaps() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))

    self.lib.fft2d_only_remaps(self.fft,cin,cout,flag)

  def only_one_remap(self,indata,outdata,flag,which):
    if "numpy" not in str(type(indata)) or "numpy" not in str(type(outdata)):
      print "Must use Numpy arrays for only_one_remap() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))
      
    self.lib.fft2d_only_one_remap(self.fft,cin,cout,flag,which)

  # tune

  def tune(self,nfast,nslow,
           in_ilo,in_ihi,in_jlo,in_jhi,
           out_ilo,out_ihi,out_jlo,out_jhi,
           permute,flag,niter,tmax,tflag):
    fftsize = c_int()
    sendsize = c_int()
    recvsize = c_int()
    self.lib.fft2d_tune(self.fft,nfast,nslow,
                        in_ilo,in_ihi,in_jlo,in_jhi,
                        out_ilo,out_ihi,out_jlo,out_jhi,
                        permute,byref(fftsize),byref(sendsize),byref(recvsize),
                        flag,niter,tmax,tflag)
    fftsize = fftsize.value
    sendsize = sendsize.value
    recvsize = recvsize.value
    return fftsize,sendsize,recvsize
    
# ------------------------------------------------------------------
# instantiate remap3dMPI thru its C-interface
# ------------------------------------------------------------------

class Remap3dMPI:

  def __init__(self,comm,precision):
    
    self.precision = precision
    if precision != 1 and precision != 2:
      print "remapMPI: precision must be 1 or 2"
      sys.exit()
    
    # load libfft3dlib.so, which contains remap2d
    
    try:
      self.lib = CDLL("libfft3dmpi.so",RTLD_GLOBAL)
    except:
      etype,value,tb = sys.exc_info()
      traceback.print_exception(etype,value,tb)
      raise OSError,"Could not load fft3dMPI dynamic library"

    # define ctypes API for each library method

    if mpi4pyflag and MPI._sizeof(comm) == sizeof(c_int): MPI_Comm = c_int
    else: MPI_Comm = c_void_p

    self.lib.remap3d_create.argtypes = [MPI_Comm,c_int,POINTER(c_void_p)]
    self.lib.remap3d_create.restype = None

    self.lib.remap3d_destroy.argtypes = [c_void_p]
    self.lib.remap3d_destroy.restype = None

    self.lib.remap3d_set.argtypes = [c_void_p,c_char_p,c_int]
    self.lib.remap3d_set.restype = None
    
    self.lib.remap3d_setup.argtypes = \
      [c_void_p,c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,c_int,c_int,c_int,c_int,c_int,
       c_int,c_int,c_int,POINTER(c_int)]
    self.lib.remap3d_setup.restype = None
    
    if precision == 1:
      self.lib.remap3d_remap.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),
           POINTER(c_float),POINTER(c_float)]
    elif precision == 2:
      self.lib.remap3d_remap.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),
           POINTER(c_double),POINTER(c_double)]
    self.lib.remap3d_compute.restype = None

    # create an instance of Remap3d

    self.remap = c_void_p()
    comm_ptr = MPI._addressof(comm)
    comm_value = MPI_Comm.from_address(comm_ptr)
    self.lib.remap3d_create(comm_value,precision,byref(self.remap))

  # destroy instance of Remap3d
  
  def __del__(self):
    self.lib.remap3d_destroy(self.remap)

  def destroy(self):
    self.lib.remap3d_destroy(self.remap)
    self.lib = None

  # setup

  def set(self,keyword,value):
    self.lib.remap3d_set(self.remap,keyword,value)

  def setup(self,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
            out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
            nqty,permute,memoryflag):
    self.lib.remap3d_setup(self.remap,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                           out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                           nqty,permute,memoryflag,sendsize,recvsize)
    return sendsize,recvsize
      
  # remap

  def remap(self,indata,outdata,sendbuf,recvbuf):
    if "numpy" not in str(type(indata)) or \
       "numpy" not in str(type(outdata)) or \
       "numpy" not in str(type(sendbuf)) or \
       "numpy" not in str(type(recvbuf)):
      print "Must use Numpy arrays for remap() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
      csendbuf = sendbuf.ctypes.data_as(POINTER(c_float))
      crecvbuf = recvbuf.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))
      csendbuf = sendbuf.ctypes.data_as(POINTER(c_double))
      crecvbuf = recvbuf.ctypes.data_as(POINTER(c_double))
      
    self.lib.remap3d_remap(self.remap,cin,cout,csendbuf,crecvbuf)

# ------------------------------------------------------------------
# instantiate remap2dMPI thru its C-interface
# ------------------------------------------------------------------

class Remap2dMPI:

  def __init__(self,comm,precision):
    
    self.precision = precision
    if precision != 1 and precision != 2:
      print "remapMPI: precision must be 1 or 2"
      sys.exit()
    
    # load libfft2dlib.so, which contains remap2d
    
    try:
      self.lib = CDLL("libfft2dmpi.so",RTLD_GLOBAL)
    except:
      etype,value,tb = sys.exc_info()
      traceback.print_exception(etype,value,tb)
      raise OSError,"Could not load fft2dMPI dynamic library"

    # define ctypes API for each library method

    if mpi4pyflag and MPI._sizeof(comm) == sizeof(c_int): MPI_Comm = c_int
    else: MPI_Comm = c_void_p

    self.lib.remap2d_create.argtypes = [MPI_Comm,c_int,POINTER(c_void_p)]
    self.lib.remap2d_create.restype = None

    self.lib.remap2d_destroy.argtypes = [c_void_p]
    self.lib.remap2d_destroy.restype = None

    self.lib.remap2d_set.argtypes = [c_void_p,c_char_p,c_int]
    self.lib.remap2d_set.restype = None
    
    self.lib.remap2d_setup.argtypes = \
    [c_void_p,c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,
        c_int,c_int,c_int,POINTER(c_int)]
    self.lib.remap2d_setup.restype = None
    
    if precision == 1:
      self.lib.remap2d_remap.argtypes = \
        [c_void_p,POINTER(c_float),POINTER(c_float),
         POINTER(c_float),POINTER(c_float)]
    elif precision == 2:
      self.lib.remap2d_remap.argtypes = \
        [c_void_p,POINTER(c_double),POINTER(c_double),
         POINTER(c_double),POINTER(c_double)]
    self.lib.remap2d_compute.restype = None

    # create an instance of Remap2d

    self.remap = c_void_p()
    comm_ptr = MPI._addressof(comm)
    comm_value = MPI_Comm.from_address(comm_ptr)
    self.lib.remap2d_create(comm_value,precision,byref(self.remap))

  # destroy instance of Remap2d
  
  def __del__(self):
    self.lib.remap2d_destroy(self.remap)

  def destroy(self):
    self.lib.remap2d_destroy(self.remap)
    self.lib = None

  # setup

  def set(self,keyword,value):
    self.lib.remap2d_set(self.remap,keyword,value)

  def setup(self,in_ilo,in_ihi,in_jlo,in_jhi,
            out_ilo,out_ihi,out_jlo,out_jhi,nqty,permute,memoryflag):
    self.lib.remap2d_setup(self.remap,in_ilo,in_ihi,in_jlo,in_jhi,
                           out_ilo,out_ihi,out_jlo,out_jhi,
                           nqty,permute,memoryflag,sendsize,recvsize)
    return sendsize,recvsize
      
  # remap

  def remap(self,indata,outdata,sendbuf,recvbuf):
    if "numpy" not in str(type(indata)) or \
       "numpy" not in str(type(outdata)) or \
       "numpy" not in str(type(sendbuf)) or \
       "numpy" not in str(type(recvbuf)):
      print "Must use Numpy arrays for remap() data"
      sys.exit()
      
    if self.precision == 1:
      cin = indata.ctypes.data_as(POINTER(c_float))
      cout = outdata.ctypes.data_as(POINTER(c_float))
      csendbuf = sendbuf.ctypes.data_as(POINTER(c_float))
      crecvbuf = recvbuf.ctypes.data_as(POINTER(c_float))
    else:
      cin = indata.ctypes.data_as(POINTER(c_double))
      cout = outdata.ctypes.data_as(POINTER(c_double))
      csendbuf = sendbuf.ctypes.data_as(POINTER(c_double))
      crecvbuf = recvbuf.ctypes.data_as(POINTER(c_double))
      
    self.lib.remap2d_remap(self.remap,cin,cout,csendbuf,crecvbuf)
