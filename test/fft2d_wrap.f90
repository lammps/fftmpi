! ISO_C_binding wrapper on fft2d C interface

module fft2d_wrap

interface
  SUBROUTINE fft2d_create(comm,PRECISION,ptr) &
          BIND(c,name='fft2d_create_fortran')
    use iso_c_binding
    integer(c_int), value :: comm
    integer(c_int), value :: precision
    type(c_ptr) :: ptr
  end subroutine fft2d_create

  subroutine fft2d_destroy(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine fft2d_destroy

  SUBROUTINE fft2d_set(ptr,keyword,value) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
    integer(c_int), value :: value
  end subroutine fft2d_set

  SUBROUTINE fft2d_get(ptr,keyword,value) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
    integer(c_int) :: value
  end subroutine fft2d_get

  SUBROUTINE fft2d_setup(ptr,nfast,nslow, &
    in_ilo,in_ihi,in_jlo,in_jhi, &
    out_ilo,out_ihi,out_jlo,out_jhi, &
    permute,fftsize,sendsize,recvsize) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    INTEGER(c_int), VALUE :: nfast,nslow
    INTEGER(c_int), VALUE :: in_ilo,in_ihi,in_jlo,in_jhi
    INTEGER(c_int), VALUE :: out_ilo,out_ihi,out_jlo,out_jhi
    INTEGER(c_int), VALUE :: permute
    INTEGER(c_int) :: fftsize,sendsize,recvsize
  end subroutine fft2d_setup

  SUBROUTINE fft2d_setup_memory_single(ptr,sendbuf,recvbuf) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    real(c_float) :: sendbuf,recvbuf
  end subroutine fft2d_setup_memory_single

  SUBROUTINE fft2d_setup_memory_double(ptr,sendbuf,recvbuf) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    real(c_double) :: sendbuf,recvbuf
  end subroutine fft2d_setup_memory_double

  SUBROUTINE fft2d_compute(ptr,in,out,flag) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), value :: in,out
    integer(c_int), value :: flag
  end subroutine fft2d_compute

  SUBROUTINE fft2d_only_1d_ffts(ptr,in,flag) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), value :: in
    integer(c_int), value :: flag
  end subroutine fft2d_only_1d_ffts

  SUBROUTINE fft2d_only_remaps(ptr,in,out,flag) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), VALUE :: in,out
    integer(c_int), value :: flag
  end subroutine fft2d_only_remaps

  SUBROUTINE fft2d_only_one_remap(ptr,in,out,flag,which) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), VALUE :: in,out
    INTEGER(c_int), VALUE :: flag,which
  end subroutine fft2d_only_one_remap

  SUBROUTINE remap2d_create(comm,ptr) &
          BIND(c,name='remap2d_create_fortran')
    use iso_c_binding
    integer(c_int), value :: comm
    type(c_ptr) :: ptr
  end subroutine remap2d_create

  subroutine remap2d_destroy(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine remap2d_destroy

  SUBROUTINE remap2d_set(ptr,keyword,value) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
    integer(c_int), value :: value
  end subroutine remap2d_set

  SUBROUTINE remap2d_setup(ptr, &
    in_ilo,in_ihi,in_jlo,in_jhi, &
    out_ilo,out_ihi,out_jlo,out_jhi, &
    nqty,permute,memoryflag,sendsize,recvsize) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    INTEGER(c_int), VALUE :: in_ilo,in_ihi,in_jlo,in_jhi
    INTEGER(c_int), VALUE :: out_ilo,out_ihi,out_jlo,out_jhi
    INTEGER(c_int), VALUE :: nqty,permute,memoryflag
    INTEGER(c_int) :: sendsize,recvsize
  end subroutine remap2d_setup

  SUBROUTINE remap2d_remap(ptr,in,out,sendbuf,recvbuf) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), VALUE :: in,out,sendbuf,recvbuf
  end subroutine remap2d_remap
end interface

end module fft2d_wrap
