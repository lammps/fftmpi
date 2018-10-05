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

! ISO_C_binding wrapper on fftMPI C interface for 3d Remaps

module remap3d_wrap

interface
  subroutine remap3d_create(comm,ptr) &
          bind(c,name='remap3d_create_fortran')
    use iso_c_binding
    integer(c_int), value :: comm
    type(c_ptr) :: ptr
  end subroutine remap3d_create

  subroutine remap3d_destroy(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine remap3d_destroy

  subroutine remap3d_set(ptr,keyword,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    character(c_char) :: keyword(*)
    integer(c_int), value :: value
  end subroutine remap3d_set

  subroutine remap3d_setup(ptr, &
    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi, &
    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi, &
    nqty,permute,memoryflag,sendsize,recvsize) &
    bind(c,name='remap3d_setup_fortran')
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
    integer(c_int), value :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
    integer(c_int), value :: nqty,permute,memoryflag
    integer(c_int) :: sendsize,recvsize
  end subroutine remap3d_setup

  subroutine remap3d_remap(ptr,in,out,sendbuf,recvbuf) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr), value :: in,out,sendbuf,recvbuf
  end subroutine remap3d_remap
end interface

end module remap3d_wrap
