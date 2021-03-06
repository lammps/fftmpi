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

! ISO_C_binding wrapper on fftMPI C interface for 3d FFTs

module fft3d_wrap

interface
  subroutine fft3d_create(comm,precision,ptr) &
          bind(c,name='fft3d_create_fortran')
    use iso_c_binding
    integer(c_int), value :: comm
    integer(c_int), value :: precision
    type(c_ptr) :: ptr
  end subroutine fft3d_create

  subroutine fft3d_destroy(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine fft3d_destroy

  subroutine fft3d_set(ptr,keyword,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
    integer(c_int), value :: value
  end subroutine fft3d_set

  function fft3d_get_int(ptr,keyword) bind(c)
    use iso_c_binding
    integer(c_int) :: fft3d_get_int
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
  end function fft3d_get_int

  function fft3d_get_double(ptr,keyword) bind(c)
    use iso_c_binding
    real(c_double) :: fft3d_get_double
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
  end function fft3d_get_double

  function fft3d_get_int64(ptr,keyword) bind(c)
    use iso_c_binding
    integer(c_int64_t) :: fft3d_get_int64
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
  end function fft3d_get_int64

  function fft3d_get_string(ptr,keyword,len) bind(c)
    use iso_c_binding
    type(c_ptr) :: fft3d_get_string
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
    integer(c_int) :: len
  end function fft3d_get_string

  function fft3d_get_int_vector(ptr,keyword,len) bind(c)
    use iso_c_binding
    type(c_ptr) :: fft3d_get_int_vector
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
    integer(c_int) :: len
  end function fft3d_get_int_vector

  function fft3d_get_double_vector(ptr,keyword,len) bind(c)
    use iso_c_binding
    type(c_ptr) :: fft3d_get_double_vector
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
    integer(c_int) :: len
  end function fft3d_get_double_vector

  subroutine fft3d_setup(ptr,nfast,nmid,nslow, &
    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi, &
    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi, &
    permute,fftsize,sendsize,recvsize) &
    bind(c,name='fft3d_setup_fortran')
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: nfast,nmid,nslow
    integer(c_int), value :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
    integer(c_int), value :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
    integer(c_int), value :: permute
    integer(c_int) :: fftsize,sendsize,recvsize
  end subroutine fft3d_setup

  subroutine fft3d_setup_memory(ptr,sendbuf,recvbuf) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr) :: sendbuf,recvbuf
  end subroutine fft3d_setup_memory

  subroutine fft3d_compute(ptr,in,out,flag) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr), value :: in,out
    integer(c_int), value :: flag
  end subroutine fft3d_compute

  subroutine fft3d_only_1d_ffts(ptr,in,flag) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr), value :: in
    integer(c_int), value :: flag
  end subroutine fft3d_only_1d_ffts

  subroutine fft3d_only_remaps(ptr,in,out,flag) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr), value :: in,out
    integer(c_int), value :: flag
  end subroutine fft3d_only_remaps

  subroutine fft3d_only_one_remap(ptr,in,out,flag,which) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr), value :: in,out
    integer(c_int), value :: flag,which
  end subroutine fft3d_only_one_remap

  subroutine fft3d_tune(ptr,nfast,nmid,nslow, &
    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi, &
    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi, &
    permute,fftsize,sendsize,recvsize,flag,niter,tmax,tflag) &
    bind(c,name='fft3d_tune_fortran')
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: nfast,nmid,nslow
    integer(c_int), value :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
    integer(c_int), value :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
    integer(c_int), value :: permute
    integer(c_int) :: fftsize,sendsize,recvsize
    integer(c_int), value :: flag,niter
    real(c_double), value :: tmax
    integer(c_int), value :: tflag
  end subroutine fft3d_tune
end interface

end module fft3d_wrap
