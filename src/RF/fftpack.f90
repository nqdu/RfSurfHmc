subroutine rfft(inp,out,n)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

    
  integer(c_int),intent(in) :: n
  !complex(c_double_complex),INTENT(IN) :: inp(n)
  real(c_double),intent(in) :: inp(n)
  complex(c_double_complex),intent(inout) :: out(n/2+1)

  ! fft plan
  type(C_PTR) :: plan
  real(c_double) ::  inp1(n)
  inp1(:) = inp(:)
  plan = fftw_plan_dft_r2c_1d(n,inp1,out,FFTW_ESTIMATE)

  call fftw_execute_dft_r2c(plan, inp1, out)
  call fftw_destroy_plan(plan)

end subroutine rfft

subroutine irfft(inp,out,n)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
    
  integer(c_int),intent(in) :: n
  complex(c_double_complex),intent(in) :: inp(n/2+1)
  real(c_double),intent(inout) :: out(n)

  ! fft plan
  type(C_PTR) :: plan
  complex(c_double_complex) inp1(n/2+1)
  inp1(:) = inp(:)
  plan = fftw_plan_dft_c2r_1d(n,inp1,out,FFTW_ESTIMATE)

  call fftw_execute_dft_c2r(plan, inp1, out)
  call fftw_destroy_plan(plan)

  out = out / n

end subroutine irfft