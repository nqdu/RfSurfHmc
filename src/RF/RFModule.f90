module RFModule
use iso_c_binding,only : c_double_complex,c_double,c_int 
implicit none

integer,parameter :: dp = c_double, dcp = c_double_complex
private
public :: cal_rf_par_time,cal_rf_time,cal_rf_freq,cal_rf_par_freq

contains 
subroutine cal_rf_par_time(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,&
                            ray_p,f0,time_shift,&
                            rf_type,par_type,rcv_fun,rcv_fun_p)&
                            bind(c,name='cal_rf_par_time_')
  implicit none
  real(dp),value,intent(in)           :: ray_p,dt,f0, time_shift
  integer(c_int), value, intent(in)   :: nlayer, rf_type, nt,par_type
  real(dp), intent(inout)             :: rcv_fun_p(nt,nlayer),rcv_fun(nt)
  real(dp), intent(in)                :: thk(nlayer), rho(nlayer),vp(nlayer), &
                                          vs(nlayer),qa(nlayer),qb(nlayer)
  
  ! local
  integer(c_int)                      :: nft, n2,  i, ilayer_inv,ilayer, it, par_layer
  complex(dcp), dimension(4,4)        :: matrix_E_inv, a_syn, a_syn_m,matrix_E_inv_par
  complex(dcp)                        :: matrix_ones(4,4),omega
  complex(dcp)                        :: imag_i,alpha(nlayer),beta(nlayer)
  complex(dcp), allocatable           :: R22(:), R21(:),R21_square(:)
  complex(dcp), allocatable           :: R22_m(:,:), R21_m(:,:),num(:)
  real(dp),parameter                  :: PI = atan(1.0) * 4.0
  real(dp),dimension(:),allocatable   :: rcv_tmp,ux,uz
  complex(dcp)                        :: all_a(4,4,nlayer),all_a_m(4,4,nlayer)

  ! padding to 2^n , total number of data point will be nft
  call nextpow2(nt,nft)
  n2 = nft / 2 + 1  

  ! allcoate space
  allocate(R22(n2), R21(n2))
  allocate( R22_m(n2,nlayer), R21_m(n2,nlayer))
  allocate(rcv_tmp(nft),ux(nft),uz(nft),R21_square(n2),num(n2))

  ! initialize identity matrix
  matrix_ones(:,:) = (0.0_dp,0.0_dp)
  do i=1,4 
    matrix_ones(i,i) = (1.0_dp,0.0_dp)
  enddo

  ! add attenuation
  imag_i = (0.0_dp,1.0_dp)
  alpha = vp * (1.0 + imag_i / (2.0 * qa) + 1.0 / (8.0 * qa**2) )
  beta = vs * (1.0 + imag_i / (2.0 * qb) + 1.0 / (8.0 * qb**2) )

 ! compute rf in frequency domain
  imag_i = (0.0_dp,1.0_dp)
  do it = 1, n2
    omega = 1.0_dp / nft / dt * (it - 1) * 2.0_dp * PI 

    ! prepare all the matirx which are needed
    do ilayer = 1, nlayer - 1
      call cal_matrix_a(omega,ray_p,thk(ilayer),alpha(ilayer),&
                        beta(ilayer),rho(ilayer),all_a(:,:,ilayer))
      call cal_matrix_a_par(omega,ray_p,thk(ilayer),alpha(ilayer),&
                        beta(ilayer),rho(ilayer),all_a_m(:,:,ilayer),par_type)
      if(par_type == 3) then
        all_a_m(:,:,ilayer) = all_a_m(:,:,ilayer) *  beta(ilayer) / vs(ilayer);
      else if(par_type == 2) then 
        all_a_m(:,:,ilayer) = all_a_m(:,:,ilayer) *  alpha(ilayer) / vp(ilayer);
      endif
    enddo

    ! prepare Einv and Einv_par
    call cal_E_inv(omega,ray_p,alpha(nlayer),beta(nlayer),&
                  rho(nlayer),matrix_E_inv)
    call cal_E_inv_par(omega,ray_p,alpha(nlayer),beta(nlayer),&
                        rho(nlayer),matrix_E_inv_par,par_type)
    if(par_type == 3) then
      matrix_E_inv_par = matrix_E_inv_par *  beta(nlayer) / vs(nlayer);
    else if(par_type == 2) then 
      matrix_E_inv_par = matrix_E_inv_par *  alpha(nlayer) / vp(nlayer);
    endif
    do par_layer = 1,nlayer 
      a_syn(:,:) = matrix_ones(:,:)
      do ilayer = 1, nlayer - 1
        ilayer_inv = nlayer - ilayer
        a_syn = matmul(a_syn,all_a(:,:,ilayer_inv))
      enddo
    enddo 
    a_syn = matmul(matrix_E_inv,a_syn)

    if (rf_type == 1) then
      R22(it) = a_syn(2,2) * imag_i
      R21(it) = a_syn(2,1) 
    else
      R22(it) = -a_syn(1,1) * imag_i
      R21(it) = a_syn(1,2) 
    endif  

    ! remove NAN
    if (isnan(abs((R22(it))))) then
      R22(it) = (0.0_dp,0.0_dp)
    endif              
    if (isnan(abs((R21(it))))) then
        R21(it) = (0.0_dp,0.0_dp)
    endif

    ! compute partial derivatives
    do par_layer = 1, nlayer
      a_syn_m(:,:) = matrix_ones(:,:)

      do ilayer = 1, nlayer - 1
        ilayer_inv = nlayer - ilayer
        if (par_layer == ilayer_inv) then
          a_syn_m = matmul(a_syn_m,all_a_m(:,:,ilayer_inv))
        else
          a_syn_m = matmul(a_syn_m,all_a(:,:,ilayer_inv))
        endif
      enddo

      if (par_layer == nlayer) then
        a_syn_m = matmul(matrix_E_inv_par,a_syn_m)
      else
        a_syn_m = matmul(matrix_E_inv,a_syn_m)
      endif

      if (rf_type == 1) then
        R22_m(it,par_layer) = a_syn_m(2,2) * imag_i  
        R21_m(it,par_layer) = a_syn_m(2,1)
      else
        R22_m(it,par_layer) = -a_syn_m(1,1) * imag_i 
        R21_m(it,par_layer) = a_syn_m(1,2)
      endif                    
      
      ! remove NAN
      if (isnan(abs(R22_m(it,par_layer)))) then
          R22_m(it,par_layer) = (0.0_dp,0.0_dp)
      endif                
      if (isnan(abs(R21_m(it,par_layer)))) then
          R21_m(it,par_layer) = (0.0_dp,0.0_dp)
      endif
    enddo ! end par_layer
  enddo ! end it

  ! compute RF 
  call irfft(R22,ux,nft)
  call irfft(R21,uz,nft)
  call deconit(ux,uz,nft,dt,time_shift,f0,rcv_tmp)
  rcv_fun(:) = rcv_tmp(1:nt)

  ! comptue R21_square
  R21_square = R21**2
  call irfft(R21_square,uz,nft)

  ! compute derivative
  do par_layer = 1, nlayer
    num = R22_m(:,par_layer) * R21(:)  -R21_m(:,par_layer) * R22(:) 
    call irfft(num,ux,nft)
    call deconit(ux,uz,nft,dt,time_shift,f0,rcv_tmp)
    rcv_fun_p(:,par_layer) = rcv_tmp(1:nt)
  enddo

  deallocate(R22,R21,R22_m, R21_m,rcv_tmp,ux,uz,R21_square,num)

end subroutine cal_rf_par_time

subroutine cal_rf_time(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,&
                        ray_p,f0,time_shift,&
                        rf_type,rcv_fun_p)bind(c,name='cal_rf_time_')
  use iso_c_binding
  implicit none
  integer(c_int),value         :: nt,rf_type,nlayer
  real(dp),value               :: dt,ray_p,f0,time_shift
  real(dp),intent(in)          :: thk(nlayer),vp(nlayer),vs(nlayer),&
                                  rho(nlayer),qa(nlayer),qb(nlayer)
  real(dp),intent(inout)       :: rcv_fun_p(nt)
  
  ! local
  complex(dcp)                 :: alpha(nlayer),beta(nlayer)
  complex(dcp), dimension(4,4) :: matrix_E_inv,a_syn, a_syn1
  integer(c_int)               :: ilayer,nft,it,n2,i,ilayer_inv
  complex(dcp),allocatable     :: R22(:),R21(:)
  real(dp),allocatable         :: ux(:),uz(:),rcv_tmp(:)
  complex(dcp)                 :: omega
  real(dp),parameter           :: pi = atan(1.0) * 4.0
  complex(dcp),parameter :: imag_i = (0.0_dp,1.0_dp)

  ! get to nextpow2 
  call nextpow2(nt,nft)
  n2 = nft / 2 + 1
  allocate(R22(n2),R21(n2),ux(nft),uz(nft),rcv_tmp(nft))

  ! add attenuation
  alpha = vp * (1.0 + imag_i / (2.0 * qa) + 1.0 / (8.0 * qa**2) )
  beta = vs * (1.0 + imag_i / (2.0 * qb) + 1.0 / (8.0 * qb**2) ) 

  ! compute R21,R22 in frequency domain
  do it = 1, n2 
    omega = (1.0/ dt / nft) * (it - 1) * 2 * pi

    ! initialize the I matrix
    a_syn = (0.0_dp,0.0_dp)
    do i = 1, 4
      a_syn(i,i) = (1.0_dp,0.0_dp)
    enddo

    ! multiple all the things from the bottom
    do ilayer = 1, nlayer - 1
      ilayer_inv = nlayer - ilayer
      call cal_matrix_a(omega,ray_p,thk(ilayer_inv),alpha(ilayer_inv),&
                        beta(ilayer_inv),rho(ilayer_inv),a_syn1)
      a_syn = matmul(a_syn,a_syn1)
    enddo

    ! compute E_inv in half space 
    call cal_E_inv(omega,ray_p,alpha(nlayer),&
                    beta(nlayer),rho(nlayer),matrix_E_inv)
    a_syn = matmul(matrix_E_inv,a_syn)

    ! compute H/Z components 
    if (rf_type == 1) then
      R22(it) = a_syn(2,2) * imag_i
      R21(it) = a_syn(2,1) 
      !print*,a_syn(2,2) - a_syn(2,3) * tmp,a_syn(2,2),a_syn(2,3) * tmp
    elseif (rf_type == 2) then
      R22(it) = -a_syn(1,1) * imag_i
      R21(it) = a_syn(1,2)
    else
      write(*,*) "error type!!"
      stop
    endif
  enddo

  ! get u/z in time domain
  call irfft(R22,ux,nft)
  call irfft(R21,uz,nft)

  ! deconvolution
  call deconit(ux,uz,nft,dt,time_shift,f0,rcv_tmp)
  rcv_fun_p(:) = rcv_tmp(1:nt)

  ! deallocate space
  deallocate(R22,R21,ux,uz,rcv_tmp)

end subroutine cal_rf_time

subroutine cal_rf_freq(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,&
                        ray_p,f0,t0,water,&
                        rf_type,rcv_fun_p)bind(c,name='cal_rf_freq_')
  use iso_c_binding
  implicit none
  integer(c_int),value         :: nt,rf_type,nlayer
  real(dp),value               :: dt,ray_p,f0,t0,water
  real(dp),intent(in)          :: thk(nlayer),vp(nlayer),vs(nlayer),&
                                  rho(nlayer),qa(nlayer),qb(nlayer)
  real(dp),intent(inout)       :: rcv_fun_p(nt)

  ! local
  complex(dcp)                 :: alpha(nlayer),beta(nlayer)
  complex(dcp), dimension(4,4) :: matrix_E_inv,a_syn, a_syn1
  integer(c_int)               :: ilayer,nft,it,n2,i,ilayer_inv
  complex(dcp),allocatable     :: R22(:),R21(:),rcv_spec(:)
  real(dp),allocatable        :: rcv_tmp(:),gauss(:),fai(:),wa(:),w(:)
  complex(dcp)                 :: omega
  real(dp)                    :: sigma,wmax
  real(dp),parameter           :: pi = atan(1.0) * 4.0
  complex(dcp),parameter :: imag_i = (0.0_dp,1.0_dp)

  ! get to nextpow2 
  call nextpow2(nt,nft)
  n2 = nft / 2 + 1
  allocate(R22(n2),R21(n2),rcv_tmp(nft),&
          gauss(n2),wa(n2),w(n2),rcv_spec(n2))

  ! add attenuation
  alpha = vp * (1.0 + imag_i / (2.0 * qa) + 1.0 / (8.0 * qa**2) )
  beta = vs * (1.0 + imag_i / (2.0 * qb) + 1.0 / (8.0 * qb**2) ) 

  ! small attenuation term
  sigma = 1.0 / (dt) / nft * 4

  ! compute R21,R22 in frequency domain
  do it = 1, n2 
    w(it) = (1.0/ dt / nft) * (it - 1) * 2 * pi
    omega = cmplx(w(it),-sigma,kind=dcp)

    ! initialize the I matrix
    a_syn = (0.0_dp,0.0_dp)
    do i = 1, 4
      a_syn(i,i) = (1.0_dp,0.0_dp)
    enddo

    ! multiple all the things from the bottom
    do ilayer = 1, nlayer - 1
      ilayer_inv = nlayer - ilayer
      call cal_matrix_a(omega,ray_p,thk(ilayer_inv),alpha(ilayer_inv),&
                        beta(ilayer_inv),rho(ilayer_inv),a_syn1)
      a_syn = matmul(a_syn,a_syn1)
    enddo

    ! compute E_inv in half space 
    call cal_E_inv(omega,ray_p,alpha(nlayer),&
                    beta(nlayer),rho(nlayer),matrix_E_inv)
    a_syn = matmul(matrix_E_inv,a_syn)

    ! compute H/Z components 
    if (rf_type == 1) then
      R22(it) = a_syn(2,2) * imag_i
      R21(it) = a_syn(2,1) 
      !print*,a_syn(2,2) - a_syn(2,3) * tmp,a_syn(2,2),a_syn(2,3) * tmp
    elseif (rf_type == 2) then
      R22(it) = -a_syn(1,1) * imag_i
      R21(it) = a_syn(1,2)
    else
      write(*,*) "error type!!"
      stop
    endif
  enddo

  ! compute gaussian 
  gauss = exp(-(w/2/f0)**2)

  ! denominator
  wa = real(R21 * conjg(R21),kind=dp)
  wmax = maxval(wa)
  fai = max(wa,water * wmax)

  ! RF in frequency domain
  rcv_spec(:) = conjg(R21) * R22 * gauss * exp(-imag_i * w * t0) / fai

  ! RF in time domain, we apply e^(sigma t ) to recover waveform
  call irfft(rcv_spec,rcv_tmp,nft)
  do it = 1,nt
    rcv_fun_p(it) = rcv_tmp(it) / dt * exp(sigma * (-t0 + (it-1) * dt)) 
  enddo

  ! free space
  deallocate(R22,R21,rcv_tmp,gauss,wa,w,rcv_spec)

end subroutine cal_rf_freq 


subroutine cal_rf_par_freq(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,&
                          ray_p,f0,t0,water,&
                          rf_type,par_type,rcv_fun,rcv_fun_p)&
                          bind(c,name='cal_rf_par_freq_')
  implicit none
  real(dp),value,intent(in)           :: ray_p,dt,f0,t0,water
  integer(c_int), value, intent(in)   :: nlayer, rf_type, nt,par_type
  real(dp), intent(inout)             :: rcv_fun_p(nt,nlayer),rcv_fun(nt)
  real(dp), intent(in)                :: thk(nlayer), rho(nlayer),vp(nlayer), &
                                          vs(nlayer),qa(nlayer),qb(nlayer)

  ! local
  integer(c_int)                      :: nft, n2,i, ilayer_inv,ilayer, it, par_layer
  complex(dcp), dimension(4,4)        :: matrix_E_inv, a_syn, a_syn_m,matrix_E_inv_par
  complex(dcp)                        :: matrix_ones(4,4)
  complex(dcp)                        :: alpha(nlayer),beta(nlayer)
  complex(dcp), allocatable           :: R22(:), R21(:),R21_square(:)
  complex(dcp), allocatable           :: R22_m(:,:), R21_m(:,:),rcv_spec(:)
  complex(dcp)                        :: all_a(4,4,nlayer),all_a_m(4,4,nlayer)
  real(dp),allocatable                :: rcv_tmp(:),gauss(:),fai(:),wa(:),w(:)
  complex(dcp)                        :: omega
  real(dp)                            :: sigma,wmax
  real(dp),parameter                  :: pi = atan(1.0) * 4.0
  complex(dcp),parameter              :: imag_i = (0.0_dp,1.0_dp)

  ! padding to 2^n , total number of data point will be nft
  call nextpow2(nt,nft)
  n2 = nft / 2 + 1  

  ! allcoate space
  allocate(R22(n2), R21(n2),rcv_tmp(nft),&
          gauss(n2),wa(n2),w(n2),rcv_spec(n2),fai(n2))
  allocate(R22_m(n2,nlayer), R21_m(n2,nlayer),R21_square(n2))

  ! initialize identity matrix
  matrix_ones(:,:) = (0.0_dp,0.0_dp)
  do i=1,4 
    matrix_ones(i,i) = (1.0_dp,0.0_dp)
  enddo

  ! add attenuation
  alpha = vp * (1.0 + imag_i / (2.0 * qa) + 1.0 / (8.0 * qa**2) )
  beta = vs * (1.0 + imag_i / (2.0 * qb) + 1.0 / (8.0 * qb**2) ) 

  ! small attenuation term
  sigma = 1.0 / (dt) / nft * 4.

  do it = 1, n2
    w(it) = 1.0_dp / nft / dt * (it - 1) * 2.0_dp * PI 
    omega = cmplx(w(it),-sigma,kind=dcp)

    ! prepare all the matirx which are needed
    do ilayer = 1, nlayer - 1
      call cal_matrix_a(omega,ray_p,thk(ilayer),alpha(ilayer),&
                        beta(ilayer),rho(ilayer),all_a(:,:,ilayer))
      call cal_matrix_a_par(omega,ray_p,thk(ilayer),alpha(ilayer),&
                        beta(ilayer),rho(ilayer),all_a_m(:,:,ilayer),par_type)
      if(par_type == 3) then
        all_a_m(:,:,ilayer) = all_a_m(:,:,ilayer) *  beta(ilayer) / vs(ilayer);
      else if(par_type == 2) then 
        all_a_m(:,:,ilayer) = all_a_m(:,:,ilayer) *  alpha(ilayer) / vp(ilayer);
      endif
    enddo

    ! prepare Einv and Einv_par
    call cal_E_inv(omega,ray_p,alpha(nlayer),beta(nlayer),&
                  rho(nlayer),matrix_E_inv)
    call cal_E_inv_par(omega,ray_p,alpha(nlayer),beta(nlayer),&
                        rho(nlayer),matrix_E_inv_par,par_type)
    if(par_type == 3) then
      matrix_E_inv_par = matrix_E_inv_par *  beta(nlayer) / vs(nlayer);
    else if(par_type == 2) then 
      matrix_E_inv_par = matrix_E_inv_par *  alpha(nlayer) / vp(nlayer);
    endif
    do par_layer = 1,nlayer 
      a_syn(:,:) = matrix_ones(:,:)
      do ilayer = 1, nlayer - 1
        ilayer_inv = nlayer - ilayer
        a_syn = matmul(a_syn,all_a(:,:,ilayer_inv))
      enddo
    enddo 
    a_syn = matmul(matrix_E_inv,a_syn)

    if (rf_type == 1) then
      R22(it) = a_syn(2,2) * imag_i
      R21(it) = a_syn(2,1) 
    else
      R22(it) = -a_syn(1,1) * imag_i
      R21(it) = a_syn(1,2) 
    endif  

    ! remove NAN
    if (isnan(abs((R22(it))))) then
      R22(it) = (0.0_dp,0.0_dp)
    endif              
    if (isnan(abs((R21(it))))) then
        R21(it) = (0.0_dp,0.0_dp)
    endif

    ! compute partial derivatives
    do par_layer = 1, nlayer
      a_syn_m(:,:) = matrix_ones(:,:)

      do ilayer = 1, nlayer - 1
        ilayer_inv = nlayer - ilayer
        if (par_layer == ilayer_inv) then
          a_syn_m = matmul(a_syn_m,all_a_m(:,:,ilayer_inv))
        else
          a_syn_m = matmul(a_syn_m,all_a(:,:,ilayer_inv))
        endif
      enddo

      if (par_layer == nlayer) then
        a_syn_m = matmul(matrix_E_inv_par,a_syn_m)
      else
        a_syn_m = matmul(matrix_E_inv,a_syn_m)
      endif

      if (rf_type == 1) then
        R22_m(it,par_layer) = a_syn_m(2,2) * imag_i  
        R21_m(it,par_layer) = a_syn_m(2,1)
      else
        R22_m(it,par_layer) = -a_syn_m(1,1) * imag_i 
        R21_m(it,par_layer) = a_syn_m(1,2)
      endif                    
      
      ! remove NAN
      if (isnan(abs(R22_m(it,par_layer)))) then
          R22_m(it,par_layer) = (0.0_dp,0.0_dp)
      endif                
      if (isnan(abs(R21_m(it,par_layer)))) then
          R21_m(it,par_layer) = (0.0_dp,0.0_dp)
      endif
    enddo ! end par_layer
  enddo ! end it

  ! compute gaussian 
  gauss = exp(-(w/2/f0)**2)

  ! denominator
  wa = real(R21 * conjg(R21),kind=dp)
  wmax = maxval(wa)
  fai = max(wa,water * wmax)

  ! RF in frequency domain
  rcv_spec(:) = conjg(R21) * R22 * gauss * exp(-imag_i * w * t0) / fai

  ! RF in time domain, we apply e^(sigma t ) to recover waveform
  call irfft(rcv_spec,rcv_tmp,nft)
  do it = 1,nt
    rcv_fun(it) = rcv_tmp(it) / dt * exp(sigma * (-t0 + (it-1) * dt)) 
  enddo

  ! comptue R21_square
  R21_square = R21**2
  wa = real(R21_square * conjg(R21_square),kind=dp)
  wmax = maxval(wa)
  fai = max(wa,water * wmax)

  ! RF derivative in frequency domain
  do par_layer = 1,nlayer
    rcv_spec = conjg(R21_square) * (R22_m(:,par_layer) * R21 - R21_m(:,par_layer) * R22) &
               *gauss * exp(-imag_i * w * t0) / fai
    call irfft(rcv_spec,rcv_tmp,nft)
    do it = 1,nt
      rcv_fun_p(it,par_layer) = rcv_tmp(it) / dt * exp(sigma * (-t0 + (it-1) * dt)) 
    enddo
  enddo

  deallocate(R22,R21,R22_m, R21_m,rcv_tmp,R21_square,w,wa,fai)
  deallocate(gauss,rcv_spec)

end subroutine cal_rf_par_freq

subroutine cal_matrix_a(omega, ray_p, thick,alpha,beta,rho,matrix_a)
  implicit none
  
  real(dp), intent(in)          :: ray_P,rho,thick 
  complex(dcp),intent(in)       :: alpha,beta,omega
  complex(dcp),intent(out)      :: matrix_a(4,4) 

  ! local
  complex(dcp)                   ::  k
  complex(dcp)                  :: miu,v_alpha, v_beta , k_alpha, k_beta, gamma, gamma1
  complex(dcp)                  :: va_k, vb_k  !~ v_alpha / k and v_beta / k 
  complex(dcp)                  :: c_a, x_a, y_a, c_b, x_b, y_b

  ! elastic modulus
  miu = rho * beta**2 
  
  ! wavenumber for vp and vs
  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2) 
  va_k = sqrt(ray_p**2 - 1.0 / alpha**2) / ray_p 
  vb_k = sqrt(ray_p**2 - 1.0  / beta**2) / ray_p

  ! gamma,gamma1
  gamma = 2 * ray_p **2 * beta ** 2
  gamma1 = 1 - 1 / gamma

  ! other locals for convenience
  c_a = cosh(v_alpha * thick)
  x_a = va_k * sinh(v_alpha * thick)
  y_a = sinh(v_alpha * thick) / va_k
  c_b = cosh(v_beta * thick)
  x_b = vb_k * sinh(v_beta * thick)
  y_b = sinh(v_beta * thick) /vb_k

  matrix_a(1,1) = c_a - gamma1 * c_b 
  matrix_a(1,2) = gamma1 * y_a - x_b
  matrix_a(1,3) = (c_b - c_a) / 2 / miu
  matrix_a(1,4) = (x_b - y_a) / 2 / miu
  matrix_a(2,1) = gamma1 * y_b - x_a
  matrix_a(2,2) = c_b - gamma1 * c_a
  matrix_a(2,3) = (x_a - y_b) / 2 / miu
  matrix_a(2,4) = (c_a - c_b) / 2 / miu
  matrix_a(3,1) = 2 * miu * gamma1 * (c_a - c_b)
  matrix_a(3,2) = 2 * miu * (gamma1**2 * y_a - x_b)
  matrix_a(3,3) = c_b - gamma1 * c_a
  matrix_a(3,4) = x_b - gamma1 * y_a
  matrix_a(4,1) = 2 * miu* (gamma1**2 * y_b - x_a) 
  matrix_a(4,2) = 2 * miu * gamma1 * (c_b - c_a)
  matrix_a(4,3) = x_a - gamma1 * y_b
  matrix_a(4,4) = c_a - gamma1 * c_b

  matrix_a = gamma * matrix_a
end subroutine cal_matrix_a

subroutine cal_matrix_a_par(omega, ray_p,thick,alpha,beta,rho,&
                                matrix_a,par_types)
  !! compute derivative of a matrix
  !! Input:
  !! omega : angular frequency
  !! par_types: =1 for rho, =2 for vp, =3 for vs, =4 for thk
  implicit none
  real(dp), intent(in)          :: ray_P,rho,thick
  complex(dcp),intent(in)       :: alpha,beta,omega 
  complex(dcp),INTENT(INOUT)    :: matrix_a(4,4)
  integer(c_int), intent(in)    :: par_types

  ! local
  complex(dcp)                  :: k,v_alpha, v_beta , k_alpha, k_beta,&
                                     gamma2, gamma3, gamma, gamma1
  complex(dcp)                  :: c_a, x_a, y_a, c_b, x_b, y_b,miu,&
                                   va_k,vb_k

  ! elastic modulus
  miu = rho * beta**2 

  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)     
  gamma = 2. * ray_p**2 * beta**2
  gamma1 = 1. - 1. / gamma
  gamma2 = gamma  / (alpha * ray_p)**2
  gamma3 = 1. / (gamma - 2)
  ! v_alpha / k and v_beta / k 
  va_k = sqrt(ray_p**2 - 1.0 / alpha**2) / ray_p 
  vb_k = sqrt(ray_p**2 - 1.0  / beta**2) / ray_p

  ! other locals
  c_a = cosh(v_alpha * thick)
  x_a = va_k * sinh(v_alpha * thick)
  y_a = sinh(v_alpha * thick) / va_k
  c_b = cosh(v_beta * thick)
  x_b = vb_k * sinh(v_beta * thick)
  y_b = sinh(v_beta * thick) /vb_k


  if (par_types == 3) then
    ! beta
    matrix_a(1,1) = 2. / beta *( gamma *(c_a - c_b)  - gamma1*k*thick*y_b ) 
    matrix_a(1,2) = 2. / beta *( gamma * (y_a - x_b) - (k * thick * c_b + y_b) )
    matrix_a(1,3) = k * thick * y_b / miu / beta
    matrix_a(1,4) = (k * thick * c_b + y_b) / miu / beta
    matrix_a(2,1) = ( (y_b - x_a) + gamma1 * gamma3 * (k*thick*c_b - y_b)  ) * 2 *gamma / beta
    matrix_a(2,2) =  2. / beta * ( gamma *(c_b - c_a) + k*thick*y_b )
    matrix_a(2,3) = -(k*thick*c_b - y_b) *gamma * gamma3 / miu / beta ! error in HSQ 2017 GJI
    matrix_a(2,4) = - matrix_a(1,3)
    matrix_a(3,1) = 4. * miu / beta *( (2*gamma - 1)*(c_a-c_b) - gamma1*k*thick*y_b )
    matrix_a(3,2) = 4. * miu / beta * ( (2*gamma)*(gamma1*y_a-x_b) - (k*thick*c_b + y_b) )
    matrix_a(3,3) =  matrix_a(2,2)
    matrix_a(3,4) = - matrix_a(1,2)
    matrix_a(4,1)= 4.*miu*gamma/beta*(2*gamma1*y_b-2*x_a + gamma1**2*gamma3*(k*thick*c_b-y_b))
    matrix_a(4,2) = - matrix_a(3,1)
    matrix_a(4,3) = - matrix_a(2,1)
    matrix_a(4,4) =  matrix_a(1,1)
  elseif(par_types == 2) then
    ! alpha
    matrix_a(1,1) = k * thick * y_a *gamma2 / alpha
    matrix_a(1,2) = 1 / va_k**2 / alpha *gamma1*gamma2 *(k*thick*c_a - y_a)
    matrix_a(1,3) = -k * thick * y_a *gamma2 / 2 / miu / alpha
    matrix_a(1,4) = - 1 / va_k**2 *(k*thick*c_a - y_a) * gamma2 / 2 / miu / alpha
    matrix_a(2,1) = - (k*thick*c_a + y_a) *gamma2 / alpha 
    matrix_a(2,2) = - k * thick * y_a * gamma1 *gamma2 / alpha
    matrix_a(2,3) = (k*thick*c_a + y_a) *gamma2 / 2 / miu / alpha
    matrix_a(2,4) = k * thick * y_a * gamma2 / 2 / miu / alpha
    matrix_a(3,1) = k * thick * y_a * gamma1 * gamma2 * 2 * miu / alpha
    matrix_a(3,2) = 2.*miu/alpha*gamma1**2 *gamma2 *(k*thick*c_a - y_a) / va_k**2
    matrix_a(3,3) = - k *thick * y_a * gamma1 *gamma2 / alpha
    matrix_a(3,4) = - 1. / alpha / va_k**2 *(k*thick*c_a - y_a) * gamma1 *gamma2
    matrix_a(4,1) = - 2. * miu / alpha *(k*thick*c_a + y_a) * gamma2
    matrix_a(4,2) = -2. * miu / alpha *k*thick*y_a*gamma1*gamma2
    matrix_a(4,3) = (k*thick*c_a + y_a) * gamma2 / alpha
    matrix_a(4,4) = k*thick*y_a / alpha
  else if(par_types == 1) then
    ! rou
    matrix_a(:,:) = (0.0_dp,0.0_dp)
    matrix_a(1,3) = - gamma / (2*rho*miu) * (-c_a + c_b) 
    matrix_a(1,4) = - gamma / (2*rho*miu) * (-y_a + x_b)
    matrix_a(2,3) = - gamma / (2*rho*miu) * (x_a - y_b)
    matrix_a(2,4) = - gamma / (2*rho*miu) * (c_a - c_b)
    matrix_a(3,1) = 2. * miu * gamma *gamma1 / rho *(c_a - c_b)
    matrix_a(3,2) = 2. * miu * gamma / rho * (gamma1**2 *y_a - x_b)
    matrix_a(4,1) = 2. * miu * gamma / rho *(-x_a + gamma1**2 * y_b)
    matrix_a(4,2) = 2. * miu * gamma * gamma1 / rho *(-c_a + c_b)
  else if(par_types == 4) then
    !thick
    matrix_a(1,1) = (x_a - gamma1*x_b) * k
    matrix_a(1,2) = gamma1*k*c_a - v_beta* vb_k *c_b
    matrix_a(1,3) = (x_b - x_a) * k / 2 / miu
    matrix_a(1,4) = (v_beta * vb_k *c_b - k*c_a) / 2 / miu
    matrix_a(2,1) = gamma1 * k * c_b - v_alpha * va_k *c_a
    matrix_a(2,2) = (x_b - gamma1*x_a) * k  
    matrix_a(2,3) = (v_alpha * va_k * c_a - k *c_b ) / 2 / miu
    matrix_a(2,4) = (x_a - x_b) * k / 2 / miu
    matrix_a(3,1) = 2. * miu * gamma1 * k *(x_a - x_b)
    matrix_a(3,2) = 2. * miu * (gamma1**2*k*c_a - v_beta * vb_k *c_b )
    matrix_a(3,3) = (x_b - gamma1 * x_a) * k
    matrix_a(3,4) = v_beta * vb_k * c_b - gamma1 * k * c_a
    matrix_a(4,1) = 2. * miu * (gamma1**2*k*c_b - v_alpha * va_k * c_a)
    matrix_a(4,2) = 2. * miu * gamma1 * k * (x_b - x_a)
    matrix_a(4,3) = v_alpha * va_k *c_a - gamma1 * k * c_b
    matrix_a(4,4) = (x_a - gamma1 * x_b) * k 
    matrix_a = gamma * matrix_a
  else 
    print*,'partypes should be one of 1,2,3,4'
    stop 
  endif        
end subroutine cal_matrix_a_par 

subroutine cal_E_inv(omega, ray_p,alpha,beta,rho,matrix_E_inv)
  implicit none
  real(dp), intent(in)               :: ray_P,rho
  complex(dcp),intent(out)           :: matrix_E_inv(4,4)
  complex(dcp),INTENT(IN)            :: omega,alpha,beta

  ! local 
  complex(dcp)                       :: k,miu,v_alpha, v_beta , k_alpha, k_beta,&
                                         gamma, gamma1,va_k,vb_k

  ! elastic modulus
  miu = rho * beta**2 
  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)     
  gamma = 2 * ray_p **2 * beta ** 2
  gamma1 = 1 - 1 / gamma
  va_k = sqrt(ray_p**2 - (1.0_dp,0.0_dp) / alpha**2) / ray_p 
  vb_k = sqrt(ray_p**2 - (1.0_dp,0.0_dp)  / beta**2) / ray_p
  
  matrix_E_inv(1,1) = -1.0
  matrix_E_inv(1,2) = -gamma1 / va_k 
  matrix_E_inv(1,3) = 1.0 / (2. * miu )
  matrix_E_inv(1,4) = 1 / (2. * miu* va_k)
  matrix_E_inv(2,1) = gamma1 / vb_k
  matrix_E_inv(2,2) = 1.0
  matrix_E_inv(2,3) = -1 / (2. * miu * vb_k)
  matrix_E_inv(2,4) = -1.0 / (2. * miu )
  matrix_E_inv(3,1) = 1.0
  matrix_E_inv(3,2) = -gamma1  / va_k
  matrix_E_inv(3,3) = -1.0 / (2. * miu)
  matrix_E_inv(3,4) = 1 / (2. * miu * va_k)
  matrix_E_inv(4,1) = -gamma1 / vb_k
  matrix_E_inv(4,2) = 1.0
  matrix_E_inv(4,3) = 1. / (2. * miu* vb_k)
  matrix_E_inv(4,4) = -1.0 / (2. * miu )

  matrix_E_inv = matrix_E_inv * 0.5 * gamma

end subroutine cal_E_inv

subroutine cal_E_inv_par(omega, ray_p,alpha,beta,rho,matrix_E_inv, par_types)
  implicit none
  real(dp), intent(in)                :: ray_P,rho
  complex(dcp),intent(in)             :: alpha,beta,omega
  complex(dcp),intent(out)            :: matrix_E_inv(4,4)
  integer(c_int), intent(in)          :: par_types  
  
  ! local
  complex(dcp)                        :: k,v_alpha, v_beta , k_alpha, k_beta,&
                                         gamma, gamma1, gamma3,miu,va_k,vb_k  

  ! elastic modulus
  miu = rho * beta**2 
  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)     
  gamma = 2 * ray_p **2 * beta ** 2
  ! v_alpha / k and v_beta / k 
  va_k = sqrt(ray_p**2 - (1.0_dp,0.0_dp) / alpha**2) / ray_p 
  vb_k = sqrt(ray_p**2 - (1.0_dp,0.0_dp)  / beta**2) / ray_p
  gamma1 = 1 - 1 / gamma
  gamma3 = 1.0_dp / (gamma - 2_dp)
  
  if (par_types == 3) then
    !beta
    matrix_E_inv(:,:) = (0.0d0,0.0d0)
    matrix_E_inv(1,1) = -1.0
    matrix_E_inv(1,2) = -1. / va_k 
    matrix_E_inv(2,1) = (1 - gamma1*gamma3) / vb_k
    matrix_E_inv(2,2) = 1.0
    matrix_E_inv(2,3) = gamma3 / 2 / miu / vb_k
    matrix_E_inv(3,1) = 1.0
    matrix_E_inv(3,2) = -1. / va_k 
    matrix_E_inv(4,1) = -matrix_E_inv(2,1)
    matrix_E_inv(4,2) = 1.0
    matrix_E_inv(4,3) = -matrix_E_inv(2,3)
    matrix_E_inv = matrix_E_inv * gamma / beta
  elseif(par_types == 1) then
    !rou
    matrix_E_inv(:,:) = (0.0d0,0.0d0)
    matrix_E_inv(1,3) = -1.0
    matrix_E_inv(1,4) = -1. / va_k 
    matrix_E_inv(2,3) = 1. / vb_k
    matrix_E_inv(2,4) = 1.0
    matrix_E_inv(3,3) = 1.0
    matrix_E_inv(3,4) = matrix_E_inv(1,4) 
    matrix_E_inv(4,3) = -matrix_E_inv(2,3)
    matrix_E_inv(4,4) = 1.0
    matrix_E_inv = matrix_E_inv * gamma / 4.0_dp / rho / miu
  elseif(par_types == 2) then
    !alpha
    matrix_E_inv(:,:) = (0.0_dp,0.0_dp)
    matrix_E_inv(1,2) = gamma1
    matrix_E_inv(1,4) = -0.5 / miu
    matrix_E_inv(3,2) = gamma1
    matrix_E_inv(3,4) = -0.5 / miu
    ! matrix_E_inv = matrix_E_inv / alpha *&
    !                k_alpha**2 * k**3 / k_beta**2 / v_alpha**3
    matrix_E_inv = matrix_E_inv * beta**2 / alpha**3 / va_k**3
  else if(par_types == 4) then 
    ! thick
    matrix_E_inv(:,:) = (0.0_dp,0.0_dp)
  else 

  endif
end
end module RFModule