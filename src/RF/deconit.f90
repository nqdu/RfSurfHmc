subroutine nextpow2(n,nout)
  implicit none
  
  integer,intent(in) :: n
  integer,intent(out) :: nout 

  nout = 1
  do while (nout < n)
    nout = nout * 2
  enddo

  return;
end subroutine nextpow2

subroutine gauss_filter(nt,dt,f0,gauss)
  implicit none
  integer,intent(in) :: nt
  real(8),intent(in) :: dt,f0 
  real(8),intent(inout) :: gauss(nt/2+1) 

  ! local
  integer :: i
  real(8),parameter :: pi = atan(1.0) * 4.
  real(8)::freq

  do i=1,nt/2+1
    freq = (i-1) /(nt * dt)
    gauss(i) = exp(-0.25 * (2 * pi * freq / f0)**2)
  enddo 

  return; 
end subroutine gauss_filter

subroutine apply_gaussian(mydata,dt,f0,nt)
   
  implicit none
  integer,intent(in) :: nt
  real(8),intent(in) :: dt,f0 
  real(8),intent(inout):: mydata(nt)

  !local 
  real(8) :: gauss(nt/2+1) 
  complex(kind=8) :: data_freq(nt/2+1)

  call gauss_filter(nt,dt,f0,gauss)
  call rfft(mydata,data_freq,nt)
  data_freq(:) = data_freq(:) * gauss 
  call irfft(data_freq,mydata,nt)

  return;

end subroutine apply_gaussian

subroutine shift_data(mydata,dt,tshift,n)
  
  implicit none
  integer,intent(in) :: n 
  real(8),intent(in) :: dt,tshift 
  real(8),intent(inout) :: mydata(n)

  ! local
  complex(kind=8) :: data_freq(n/2+1)
  integer :: i 
  real(8),parameter :: pi = atan(1.0) * 4.0

  call rfft(mydata,data_freq,n)
  do i=1,n/2+1
    data_freq(i) = data_freq(i) * exp(-cmplx(0,1.0) * (i-1) / (n*dt) * pi * 2 * tshift)
  enddo
  call irfft(data_freq,mydata,n)

end subroutine shift_data

subroutine shift_data1(mydata,dt,tshift,n)
  
  implicit none
  integer,intent(in) :: n 
  real(8),intent(in) :: dt,tshift 
  real(8),intent(inout) :: mydata(n)

  ! local
  real(8) :: data_pad(n*2)
  complex(kind=8) :: data_freq(n+1)
  integer :: i,n2 
  real(8),parameter :: pi = atan(1.0) * 4.0
  n2 = n * 2

  data_pad = 0.
  data_pad(1:n) = mydata(:)
  call rfft(data_pad,data_freq,n2)
  do i=1,n+1
    data_freq(i) = data_freq(i) * exp(-cmplx(0,1.0) * (i-1) / (n2*dt) * pi * 2 * tshift)
  enddo
  call irfft(data_freq,data_pad,n2)
  mydata = data_pad(1:n)

end subroutine shift_data1

subroutine mycorrelate(a,b,out,n)
  
  implicit none
  
  integer,intent(in) :: n 
  real(8),intent(in) :: a(n),b(n)
  real(8),intent(inout) :: out(n)

  !local 
  complex(kind=8) :: aft(n/2+1),bft(n/2+1),c(n/2+1)
  call rfft(a,aft,n)
  call rfft(b,bft,n)
  c = aft * conjg(bft)
  
  call irfft(c,out,n)

end subroutine mycorrelate

subroutine myconvolve(a,b,out,n)
  
  implicit none
  
  integer,intent(in) :: n 
  real(8),intent(in) :: a(n),b(n)
  real(8),intent(inout) :: out(n)

  !local 
  complex(kind=8) :: aft(n/2+1),bft(n/2+1),c(n/2+1)
  call rfft(a,aft,n)
  call rfft(b,bft,n)
  c = aft * bft
  
  call irfft(c,out,n)

end subroutine myconvolve

subroutine deconit(u, w, nt,dt, tshift, f0,out)
  implicit none
  integer,intent(in) :: nt 
  real(8),intent(in) :: dt,tshift,f0,u(nt),w(nt)
  real(8),intent(inout) :: out(nt)

  ! local
  integer :: nft 
  real(8),dimension(:),allocatable :: uflt,wflt,wcopy,p,rflt 
  real(8),dimension(:),allocatable :: cuw,temp1,temp2 
  real(8) :: invpw,invpu,minderr,d_error,sumsq_i,sumsq
  integer,parameter :: maxiter = 200
  integer :: it,idx(1)

  ! allocate space 
  call nextpow2(nt,nft)
  allocate(uflt(nft),wflt(nft),wcopy(nft),p(nft),rflt(nft))
  allocate(cuw(nft),temp1(nft),temp2(nft))

  ! copy 
  wflt(:) = 0.; uflt(:) = 0.
  wflt(1:nt) = w(:); uflt(1:nt) = u(:)
  wcopy = wflt

  ! filter input arrays
  call apply_gaussian(uflt,dt,f0,nft)
  call apply_gaussian(wflt,dt,f0,nft)

  ! init 
  invpw = 1. / sum(wflt**2) / dt 
  invpu = 1. / sum(uflt**2) / dt 
  p(:) = 0.
  sumsq_i = 1.0
  sumsq = 50
  minderr = 0.001
  d_error = 100 * invpw + minderr
  rflt(:) = uflt(:)

  ! iterative deconvolution
  do it = 1,maxiter
    if( abs(d_error) <= minderr) exit 
    call mycorrelate(rflt,wflt,cuw,nft)
    cuw(:) = cuw(:) * dt
    idx = maxloc(abs(cuw(1:nft/2)))
    P(idx(1)) = P(idx(1)) + cuw(idx(1)) * invpw / dt 
    temp1 = P 
    call apply_gaussian(temp1,dt,f0,nft)
    call myconvolve(temp1,wcopy,temp2,nft)
    rflt = uflt - temp2 * dt 

    ! compute error 
    sumsq = sum(rflt**2) * dt * invpu
    d_error = 100. * (sumsq_i - sumsq)
    sumsq_i = sumsq
  enddo

  ! get rf and time shift 
  call apply_gaussian(p,dt,f0,nft)
  call shift_data(p,dt,tshift,nft)
  out(1:nt) = p(1:nt)

  deallocate(uflt,wflt,wcopy,p,rflt, cuw,temp1,temp2 )

end subroutine deconit