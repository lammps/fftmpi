! ISO_C_binding wrapper on fft3d C interface

module fft3d_wrap

interface
  SUBROUTINE fft3d_create(comm,PRECISION,ptr) &
          BIND(c,name='fft3d_create_fortran')
    use iso_c_binding
    integer(c_int), value :: comm
    integer(c_int), value :: precision
    type(c_ptr) :: ptr
  end subroutine fft3d_create

  subroutine fft3d_destroy(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine fft3d_destroy

  SUBROUTINE fft3d_set(ptr,keyword,value) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
    integer(c_int), value :: value
  end subroutine fft3d_set

  function fft3d_get_int(ptr,keyword) BIND(c)
    use iso_c_binding
    integer(c_int) :: fft3d_get_int
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
  end function fft3d_get_int

  function fft3d_get_double(ptr,keyword) BIND(c)
    use iso_c_binding
    real(c_double) :: fft3d_get_double
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
  end function fft3d_get_double

  function fft3d_get_int64(ptr,keyword) BIND(c)
    use iso_c_binding
    integer(c_int64_t) :: fft3d_get_int64
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
  end function fft3d_get_int64

  function fft3d_get_string(ptr,keyword) BIND(c)
    use iso_c_binding
    type(c_ptr) :: fft3d_get_string
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
  end function fft3d_get_string

  function fft3d_get_int_vector(ptr,keyword) BIND(c)
    use iso_c_binding
    type(c_ptr) :: fft3d_get_int_vector
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
  end function fft3d_get_int_vector

  function fft3d_get_double_vector(ptr,keyword) BIND(c)
    use iso_c_binding
    type(c_ptr) :: fft3d_get_double_vector
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
  end function fft3d_get_double_vector

  SUBROUTINE fft3d_setup(ptr,nfast,nmid,nslow, &
    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi, &
    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi, &
    permute,fftsize,sendsize,recvsize) &
    BIND(c,name='fft3d_setup_fortran')
    use iso_c_binding
    type(c_ptr), value :: ptr
    INTEGER(c_int), VALUE :: nfast,nmid,nslow
    INTEGER(c_int), VALUE :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
    INTEGER(c_int), VALUE :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
    INTEGER(c_int), VALUE :: permute
    INTEGER(c_int) :: fftsize,sendsize,recvsize
  end subroutine fft3d_setup

  SUBROUTINE fft3d_setup_memory(ptr,sendbuf,recvbuf) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr) :: sendbuf,recvbuf
  end subroutine fft3d_setup_memory

  SUBROUTINE fft3d_compute(ptr,in,out,flag) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    type(c_ptr), value :: in,out
    integer(c_int), value :: flag
  end subroutine fft3d_compute

  SUBROUTINE fft3d_only_1d_ffts(ptr,in,flag) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), value :: in
    integer(c_int), value :: flag
  end subroutine fft3d_only_1d_ffts

  SUBROUTINE fft3d_only_remaps(ptr,in,out,flag) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), VALUE :: in,out
    integer(c_int), value :: flag
  end subroutine fft3d_only_remaps

  SUBROUTINE fft3d_only_one_remap(ptr,in,out,flag,which) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), VALUE :: in,out
    INTEGER(c_int), VALUE :: flag,which
  end subroutine fft3d_only_one_remap

  SUBROUTINE fft3d_tune(ptr,nfast,nmid,nslow, &
    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi, &
    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi, &
    permute,fftsize,sendsize,recvsize,flag,niter,tmax,tflag) &
    BIND(c,name='fft3d_tune_fortran')
    use iso_c_binding
    type(c_ptr), value :: ptr
    INTEGER(c_int), VALUE :: nfast,nmid,nslow
    INTEGER(c_int), VALUE :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
    INTEGER(c_int), VALUE :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
    INTEGER(c_int), VALUE :: permute
    INTEGER(c_int) :: fftsize,sendsize,recvsize
    INTEGER(c_int), VALUE :: flag,niter
    REAL(c_double), value :: tmax
    INTEGER(c_int), VALUE :: tflag
  end subroutine fft3d_tune

  SUBROUTINE remap3d_create(comm,ptr) &
          BIND(c,name='remap3d_create_fortran')
    use iso_c_binding
    integer(c_int), value :: comm
    type(c_ptr) :: ptr
  end subroutine remap3d_create

  subroutine remap3d_destroy(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine remap3d_destroy

  SUBROUTINE remap3d_set(ptr,keyword,value) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    CHARACTER(c_char) :: keyword(*)
    integer(c_int), value :: value
  end subroutine remap3d_set

  SUBROUTINE remap3d_setup(ptr, &
    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi, &
    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi, &
    nqty,permute,memoryflag,sendsize,recvsize) &
    BIND(c,name='remap3d_setup_fortran')
    use iso_c_binding
    type(c_ptr), value :: ptr
    INTEGER(c_int), VALUE :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
    INTEGER(c_int), VALUE :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
    INTEGER(c_int), VALUE :: nqty,permute,memoryflag
    INTEGER(c_int) :: sendsize,recvsize
  end subroutine remap3d_setup

  SUBROUTINE remap3d_remap(ptr,in,out,sendbuf,recvbuf) BIND(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    TYPE(c_ptr), VALUE :: in,out,sendbuf,recvbuf
  end subroutine remap3d_remap
end interface

end module fft3d_wrap
