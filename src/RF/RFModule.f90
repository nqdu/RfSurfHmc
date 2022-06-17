module layer_var
use, intrinsic :: iso_c_binding
implicit none
integer(c_int),PARAMETER     :: sp=c_float,dp=c_double 
integer(c_int),PARAMETER     :: scp=c_float_complex,dcp=c_double_complex                        
real(dp), parameter          :: PI = atan(1.0_dp) * 4.0_dp

! GLOBAL Variables
real(dp), allocatable        :: thk(:),alpha(:),beta(:),rho(:) 

contains 
! constructor
subroutine layer_var_alloc(h,a,b,r,n)
  implicit none
  integer(c_int)        :: n 
  real(dp)              :: h(n),a(n),b(n),r(n)

  ! allocate model
  ALLOCATE(thk(n),alpha(n),beta(n),rho(n))
  thk = h; alpha = a;
  beta = b; rho = r;
end subroutine layer_var_alloc

! destructor
subroutine layer_var_dealloc()
  implicit none
  DEALLOCATE(thk,alpha,beta,rho)
end subroutine layer_var_dealloc

end module layer_var ! end module

!! Receiver Function Module
module RFModule

use,intrinsic             :: iso_c_binding
use layer_var,only        : sp,dp,dcp,scp,pi
implicit none

contains 
subroutine receiverFunc_par(thk,vp,vs,rho,nlayer,nt,dt,&
                            ray_p,gauss,time_shift,water_level,&
                            rf_type,par_type,rcvfun,rcv_fun_p)&
                            bind(c,name='receiverFunc_par_')
  use layer_var,only                   : layer_var_alloc,layer_var_dealloc
  implicit none
  real(dp),value,intent(in)           :: ray_p,dt, gauss, time_shift, water_level
  integer(c_int), value, intent(in)   :: nlayer, rf_type, nt,par_type
  real(dp), intent(inout)             :: rcv_fun_p(nt,nlayer),rcvfun(nt)
  real(dp), intent(in)                :: thk(nlayer), rho(nlayer),vp(nlayer), vs(nlayer)

  ! set variables
  call layer_var_alloc(thk,vp,vs,rho,nlayer)

  ! compute rf and its derivative
  call rf_cal_partial(rf_type,nlayer,ray_p,nt,dt,gauss,time_shift,&
                      water_level,par_type,rcvfun,rcv_fun_p)

  ! free space
  call layer_var_dealloc()
end subroutine receiverFunc_par

subroutine receiverFunc(thk,vp,vs,rho,nlayer,nt,dt,&
                        ray_p,gauss,time_shift,water_level,&
                        rf_type,rcv_fun_p)bind(c,name='receiverFunc_')
  !bind(c,name="receiverFunc_cal_")
  use layer_var,only                   : layer_var_alloc,layer_var_dealloc
  implicit none
  real(dp), value, intent(in)          :: ray_p, dt, gauss, time_shift, water_level
  integer(c_int), value, intent(in)    :: nlayer, rf_type,nt
  real(dp), intent(inout)              :: rcv_fun_p(nt)
  real(dp), intent(in)                 :: thk(nlayer), rho(nlayer), vp(nlayer), vs(nlayer)

  ! set variables
  call layer_var_alloc(thk,vp,vs,rho,nlayer)
  
  ! rf
  call rf_cal(rf_type,nlayer,ray_p,nt,dt,gauss,time_shift,&
              water_level,rcv_fun_p)

  ! free space
  call layer_var_dealloc()       
end subroutine receiverFunc

subroutine rf_cal(rf_type, nlayer, ray_p, nt, dt,&
                 gauss, time_shift,water_level, rcv_fun )
! rf_call is used for calculating  receiver function 
! the in/out parameters are as below
! Input parameter:
!   rf_type     :   (int)       type of receiver function; 1 = p wave 2 = s wave
!   nlayer      :   (int)       how many layers are calculated
!   ray_p       :   (c_double)  ray parameter in s/km, which is used for calculate the incident angel
!   nt:         :   (int)       the length of rf (number)
!   dt          :   (c_double)  the sample interval of rf (second)
!   gauss       :   (c_double)  gauss factor for filter, usually 1.5(~0.5Hz) or 3(~1Hz)
!   water_level :   (c_double)  avoiding unstable for division in frequency domain, usually choose 1e-3
! Output parameter:
!   rcv_fun  :   (c_double)  the result of rf derivative, size:(ceiling(duation/dt))
  use layer_var,only              : thk,alpha,beta,rho             
  implicit none
  real(dp), intent(in)           :: ray_p, dt, gauss, time_shift, water_level
  integer(c_int), intent(in)     :: nlayer,nt
  real(dp),INTENT(INOUT)         :: rcv_fun(nt)
  
  !local
  integer(c_int)                 :: nft, i, ilayer_inv, n2, ilayer, it, rf_type
  complex(dcp), dimension(4,4)   :: matrix_E_inv,a_syn, a_syn1
  real(dp), allocatable          :: omega(:), fai(:), R21_abs(:),rcv_fun_p_inv(:)
  complex(dcp)                   :: imag_i,cs
  complex(dcp), allocatable      :: rcv_fun_p_spectrum(:), R22(:), R21(:)
  real(dp)                       :: max_r21


  imag_i = (0.0_dp,1.0_dp)
  ! padding to 2^n , total number of data point will be nft
  nft = nextpow2(nt)
  n2 = nft / 2 + 1
  
  allocate( rcv_fun_p_spectrum(n2),R22(n2), R21(n2), fai(n2), R21_abs(n2) )
  allocate( omega(n2),rcv_fun_p_inv(nft))

  ! only perform half spectrum
  do it = 1, n2
    omega(it) = (1.0/ dt / nft) * (it - 1) * 2 * PI

    ! initialize the I matrix
    a_syn = (0.0_dp,0.0_dp)
    do i = 1, 4
      a_syn(i,i) = (1.0_dp,0.0_dp)
    enddo

    ! multiple all the things from the bottom
    do ilayer = 1, nlayer - 1
      ilayer_inv = nlayer - ilayer
      call cal_matrix_a(omega(it),ray_p,thk(ilayer_inv),alpha(ilayer_inv),&
                        beta(ilayer_inv),rho(ilayer_inv),a_syn1)
      a_syn = matmul(a_syn,a_syn1)
    enddo

    ! compute E_inv in half space 
    call cal_E_inv(omega(it),ray_p,alpha(nlayer),&
                    beta(nlayer),rho(nlayer),matrix_E_inv)
    a_syn = matmul(matrix_E_inv,a_syn)

    ! compute H/Z components 
    if (rf_type == 1) then
      R22(it) = a_syn(2,2) 
      R21(it) = a_syn(2,1) 
    elseif (rf_type == 2) then
      R22(it) = -a_syn(1,1) 
      R21(it) = a_syn(1,2)
    else
      write(*,*) "error type!!"
      stop
    endif

    !check NAN
    if (isnan(abs(R22(it)))) then
        R22(it) = (0.0_dp,0.0_dp)
    endif
    if (isnan(abs(R21(it)))) then
        R21(it) = (0.0_dp,0.0_dp)
    endif
    ! remove all the nan point
    ! inspired by nanqiao.
    ! but Junliu thinks that it is used to deal with zero freq
    ! we dont need to do that
  enddo

  R21_abs = real( R21* conjg(R21) )
  max_r21 = maxval(R21_abs)
  fai = max(R21_abs,water_level * max_r21)
  ! That's how water level used

  !compute rf in frequency domain
  do it = 1, n2
    ! cs is for time shift
    if (rf_type == 1) then
      cs= - imag_i*omega(it)*time_shift
      rcv_fun_p_spectrum(it) = imag_i*(R22(it)*conjg(R21(it)))/fai(it)&
                               *exp(-omega(it)**2/(4*gauss**2))*exp(cs)
    else
      !cs= - imag_i*omega(it)*(duration - time_shift)
      ! change the cs when calculate the srf 
      cs= imag_i*omega(it)*time_shift

      ! contrasted, we multiple the conjg 
      ! usually we do not perfrom the complex divide
      rcv_fun_p_spectrum(it) = imag_i*conjg(R22(it)*conjg(R21(it)))/fai(it)&
                                *exp(-omega(it)**2/(4*gauss**2))*exp(cs)
  
      ! there are something wrongs in s_rf calculation
      ! according to lupei zhu's program
      ! cngtv(conjg(cmltp(cmltp(IMAGE,dis.ps), cinvs(dis.ss))));
      ! it seems that the conjg is needed               
    endif
  enddo

  ! inverse ifft
  ! remember the definition of fft, check the operator of forward fft
  ! is that e^(-iwt)? 
  ! do we need to divide 1/n when we perform the invert fft? 
  call irfft(rcv_fun_p_spectrum,rcv_fun_p_inv,nft)   

  ! output result
  rcv_fun(1:nt) = rcv_fun_p_inv(1:nt) / dt

  ! free space
  DEALLOCATE(omega,fai,R21,R22,R21_abs,rcv_fun_p_inv,rcv_fun_p_spectrum)

end subroutine

subroutine rf_cal_partial(rf_type, nlayer,ray_p, nt, dt, gauss, &
                          time_shift, water_level, par_type,&
                          rcvfun,rcv_fun_p)
  !! rf_cal_partial is used for calculating (analytic partial derivative of) receiver function 
  !! the in/out parameters are as below
  !! Input parameter:
  !!   rf_type     :   (int)       type of receiver function; 1 = p wave 2 = s wave
  !!   nlayer      :   (int)       how many layers are calculated
  !!   ray_p       :   (c_double)  ray parameter in s/km, which is used for calculate the incident angel
  !!   nt:         :   (int)       the length of rf (number)
  !!   dt          :   (c_double)  the sample interval of rf (second)
  !!   gauss       :   (c_double)  gauss factor for filter, usually 1.5(~0.5Hz) or 3(~1Hz)
  !!   water_level :   (c_double)  avoiding unstable for division in frequency domain, usually choose 1e-3
  !!   par_type    :   (int)       which layer should be calculate for analytic partial derivative
  !! Output parameter:
  !!   rcvfun      :    (c_double)  receiver function
  !!   rcv_fun_p   :   (c_double)  the result of rf derivative, size:(nt,nlayer)
  use layer_var,only                   : thk,beta,alpha,rho
  implicit none
  real(dp),intent(in)                 :: ray_p, dt, gauss, time_shift, water_level
  integer(c_int), intent(in)          :: nlayer, rf_type,nt,par_type
  real(dp),intent(inout)              :: rcv_fun_p(nt,nlayer),rcvfun(nt)

  ! local 
  integer(c_int)                      :: nft, n2,  i, ilayer_inv,ilayer, it, par_layer
  complex(dcp), dimension(4,4)        :: matrix_E_inv, a_syn, a_syn_m,matrix_E_inv_par
  real(dp), allocatable               :: omega(:), fai(:), R21_square_abs(:)
  complex(dcp)                        :: matrix_ones(4,4)
  complex(dcp)                        :: imag_i, numer,cs
  complex(dcp), allocatable           :: rcv_fun_p_spectrum(:),R22(:), R21(:)
  complex(dcp), allocatable           :: R22_m(:,:), R21_m(:,:), R21_square(:)
  real(dp)                            :: max_r21
  real(dp),allocatable                :: rcv_fun_p_inv(:)
  complex(dcp)                        :: all_a(4,4,nlayer),all_a_m(4,4,nlayer)

  imag_i = (0.0_dp,1.0_dp)

  ! padding to 2^n , total number of data point will be nft
  nft = nextpow2(nt)
  n2 = nft / 2 + 1

  ! allocate space   
  allocate( rcv_fun_p_spectrum(n2),R22(n2), R21(n2), fai(n2) )
  allocate( R22_m(n2,nlayer), R21_m(n2,nlayer), R21_square(n2), R21_square_abs(n2) )
  allocate( omega(n2), rcv_fun_p_inv(nft))

  ! prepare the I matrix
  matrix_ones(:,:) = (0.0_dp,0.0_dp)
  do i=1,4 
    matrix_ones(i,i) = (1.0_dp,0.0_dp)
  enddo
      
  ! compute R21/R22, R21_m/R22_m
  do it = 1, n2
    omega(it) = 1.0_dp / nft / dt * (it - 1) * 2.0_dp * PI 

    ! prepare all the matirx which are needed
    do ilayer = 1, nlayer - 1
      call cal_matrix_a(omega(it),ray_p,thk(ilayer),alpha(ilayer),&
                        beta(ilayer),rho(ilayer),all_a(:,:,ilayer))
      call cal_matrix_a_par(omega(it),ray_p,thk(ilayer),alpha(ilayer),&
                        beta(ilayer),rho(ilayer),all_a_m(:,:,ilayer),par_type)
    enddo

    ! prepare Einv and Einv_par
    call cal_E_inv(omega(it),ray_p,alpha(nlayer),beta(nlayer),&
                  rho(nlayer),matrix_E_inv)
    call cal_E_inv_par(omega(it),ray_p,alpha(nlayer),beta(nlayer),&
                        rho(nlayer),matrix_E_inv_par,par_type)

    do par_layer = 1,nlayer 
      a_syn(:,:) = matrix_ones(:,:)
      do ilayer = 1, nlayer - 1
        ilayer_inv = nlayer - ilayer
        a_syn = matmul(a_syn,all_a(:,:,ilayer_inv))
      enddo
    enddo 
    a_syn = matmul(matrix_E_inv,a_syn)

    if (rf_type == 1) then
      R22(it) = a_syn(2,2)
      R21(it) = a_syn(2,1) 
    else
      R22(it) = -a_syn(1,1)
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
        R22_m(it,par_layer) = a_syn_m(2,2)  
        R21_m(it,par_layer) = a_syn_m(2,1)
      else
        R22_m(it,par_layer) = -a_syn_m(1,1) 
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

  
  ! water-level regularization  for rf
  max_r21 = maxval(real( R21* conjg(R21)))
  fai = max(real( R21* conjg(R21) ),water_level * max_r21)

  ! compute rf
  do it=1,n2
    !compute rf in frequency domain
    ! cs is for time shift
    if (rf_type == 1) then
      cs = - imag_i*omega(it)*time_shift
      rcv_fun_p_spectrum(it) = imag_i*(R22(it)*conjg(R21(it)))/fai(it)&
                                *exp(-omega(it)**2/(4*gauss**2))*exp(cs)
    else
      !cs= - imag_i*omega(it)*(duration - time_shift)
      ! change the cs when calculate the srf 
      cs= imag_i*omega(it)*time_shift
      rcv_fun_p_spectrum(it) = imag_i*conjg(R22(it)*conjg(R21(it)))/fai(it)&
                                *exp(-omega(it)**2/(4*gauss**2))*exp(cs)          
    endif
  enddo
  call irfft(rcv_fun_p_spectrum,rcv_fun_p_inv,nft)
  rcvfun(:) = rcv_fun_p_inv(1:nt) / dt

  ! water-level regularization for derivative
  R21_square = R21(:) **2
  R21_square_abs = real( R21_square* conjg(R21_square) )
  max_r21 = maxval(R21_square_abs)
  fai = max(R21_square_abs,water_level * max_r21)

  ! compute derivative
  do par_layer = 1, nlayer
    do it = 1, n2
      if (rf_type == 1) then
        cs = cdexp(- imag_i*omega(it)*time_shift)
        numer = imag_i*(R22_m(it,par_layer)*R21(it)- &
                R21_m(it,par_layer)*R22(it))*conjg(R21_square(it))
      else 
        cs = cdexp(imag_i*omega(it)*time_shift)
        numer = imag_i*conjg((R22_m(it,par_layer)*R21(it)- &
                R21_m(it,par_layer)*R22(it))*conjg(R21_square(it)))
      endif                    
      rcv_fun_p_spectrum(it) = numer / fai(it) * exp( -omega(it)**2/(4*gauss**2))*cs
    enddo

    ! ifft 
    call irfft(rcv_fun_p_spectrum,rcv_fun_p_inv,nft)
    rcv_fun_p(:,par_layer) = rcv_fun_p_inv(1:nt) / dt
  enddo

  DEALLOCATE(rcv_fun_p_spectrum,R22, R21, fai)
  deallocate( R22_m, R21_m, R21_square, R21_square_abs)
  deallocate( omega, rcv_fun_p_inv)
end subroutine rf_cal_partial


function nextpow2(nt) result(nft)
  implicit none
  integer(c_int),INTENT(IN)          :: nt 
  INTEGER(c_int)                     :: nft 

  nft = 1
  do while (nft < nt)
      nft = nft * 2
  end do

  return;

end function nextpow2

subroutine cal_matrix_a(omega, ray_p, thick,alpha,beta,rho,matrix_a)
  implicit none
  
  real(dp), intent(in)          :: ray_P, omega, alpha,beta,rho,thick 
  complex(dcp),intent(out)      :: matrix_a(4,4) 

  ! local
  real(dp)                      ::  k,lambda,miu
  complex(dcp)                  :: v_alpha, v_beta , k_alpha, k_beta, gamma, gamma1
  complex(dcp)                  :: c_a, x_a, y_a, c_b, x_b, y_b

  ! elastic modulus
  miu = rho * beta**2 
  lambda = rho * alpha**2 - 2. * miu
  
  ! wavenumber for vp and vs
  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)    

  ! gamma,gamma1
  gamma = 2 * k**2 * beta**2 / omega**2
  gamma1 = 1 - 1 / gamma

  ! other locals for convenience
  c_a = cosh(v_alpha * thick)
  x_a = v_alpha * sinh(v_alpha * thick) / k
  y_a = k * sinh(v_alpha * thick) / v_alpha
  c_b = cosh(v_beta * thick)
  x_b = v_beta * sinh(v_beta * thick) / k
  y_b = k * sinh(v_beta * thick) /v_beta

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
  real(dp), intent(in)          :: ray_P, omega, alpha,beta,rho,thick  
  complex(dcp),INTENT(INOUT)    :: matrix_a(4,4)
  integer(c_int), intent(in)    :: par_types

  ! local
  real(dp)                      :: k,lambda,miu
  complex(dcp)                  :: v_alpha, v_beta , k_alpha, k_beta,&
                                     gamma2, gamma3, gamma, gamma1
  complex(dcp)                  :: c_a, x_a, y_a, c_b, x_b, y_b

  ! elastic modulus
  miu = rho * beta**2 
  lambda = rho * alpha**2 - 2. * miu

  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)     
  gamma = 2. * k**2 * beta**2 / omega**2
  gamma1 = 1. - 1. / gamma
  gamma2 = gamma * k_alpha**2 / k**2
  gamma3 = 1. / (gamma - 2)

  ! other locals
  c_a = cosh(v_alpha * thick)
  x_a = v_alpha * sinh(v_alpha * thick) / k
  y_a = k * sinh(v_alpha * thick) / v_alpha
  c_b = cosh(v_beta * thick)
  x_b = v_beta * sinh(v_beta * thick) / k
  y_b = k * sinh(v_beta * thick) /v_beta

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
    matrix_a(1,2) = k**2 / v_alpha**2 / alpha *gamma1*gamma2 *(k*thick*c_a - y_a)
    matrix_a(1,3) = -k * thick * y_a *gamma2 / 2 / miu / alpha
    matrix_a(1,4) = - k**2 / v_alpha**2 *(k*thick*c_a - y_a) * gamma2 / 2 / miu / alpha
    matrix_a(2,1) = - (k*thick*c_a + y_a) *gamma2 / alpha 
    matrix_a(2,2) = - k * thick * y_a * gamma1 *gamma2 / alpha
    matrix_a(2,3) = (k*thick*c_a + y_a) *gamma2 / 2 / miu / alpha
    matrix_a(2,4) = k * thick * y_a * gamma2 / 2 / miu / alpha
    matrix_a(3,1) = k * thick * y_a * gamma1 * gamma2 * 2 * miu / alpha
    matrix_a(3,2) = 2.*miu/alpha*gamma1**2 *gamma2 *(k*thick*c_a - y_a) *k**2 / v_alpha**2
    matrix_a(3,3) = - k *thick * y_a * gamma1 *gamma2 / alpha
    matrix_a(3,4) = - 1. / alpha *k**2 / v_alpha**2 *(k*thick*c_a - y_a) * gamma1 *gamma2
    matrix_a(4,1) = - 2. * miu / alpha *(k*thick*c_a + y_a) * gamma2
    matrix_a(4,2) = -2. * miu / alpha *k*thick*y_a*gamma1*gamma2
    matrix_a(4,3) = (k*thick*c_a + y_a) * gamma2 / alpha
    matrix_a(4,4) = k*thick*y_a / alpha * gamma2 
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
    matrix_a(1,2) = gamma1*k*c_a - v_beta**2*c_b/k
    matrix_a(1,3) = (x_b - x_a) * k / 2 / miu
    matrix_a(1,4) = (v_beta**2*c_b/k - k*c_a) / 2 / miu
    matrix_a(2,1) = gamma1 * k * c_b - v_alpha**2 *c_a / k
    matrix_a(2,2) = (x_b - gamma1*x_a) * k  
    matrix_a(2,3) = (v_alpha**2 * c_a/k - k *c_b ) / 2 / miu
    matrix_a(2,4) = (x_a - x_b) * k / 2 / miu
    matrix_a(3,1) = 2. * miu * gamma1 * k *(x_a - x_b)
    matrix_a(3,2) = 2. * miu * (gamma1**2*k*c_a - v_beta**2 *c_b / k)
    matrix_a(3,3) = (x_b - gamma1 * x_a) * k
    matrix_a(3,4) = v_beta **2 * c_b / k - gamma1 * k * c_a
    matrix_a(4,1) = 2. * miu * (gamma1**2*k*c_b - v_alpha**2*c_a/k)
    matrix_a(4,2) = 2. * miu * gamma1 * k * (x_b - x_a)
    matrix_a(4,3) = v_alpha**2*c_a / k - gamma1 * k * c_b
    matrix_a(4,4) = (x_a - gamma1 * x_b) * k 
    matrix_a = gamma * matrix_a
  else 
    print*,'partypes should be one of 1,2,3,4'
    stop 
  endif        
end subroutine cal_matrix_a_par 

subroutine cal_E_inv(omega, ray_p,alpha,beta,rho,matrix_E_inv)
  implicit none
  real(dp), intent(in)               :: ray_P, omega
  complex(dcp),intent(out)           :: matrix_E_inv(4,4)
  real(dp),INTENT(IN)                :: alpha,beta,rho

  ! local
  real(dp)                           ::  k,miu   
  complex(dcp)                       :: v_alpha, v_beta , k_alpha, k_beta,&
                                         gamma, gamma1

  ! elastic modulus
  miu = rho * beta**2 
  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)     
  gamma = 2 * k**2 * beta**2 / omega**2
  gamma1 = 1 - 1 / gamma

  
  matrix_E_inv(1,1) = -1.0
  matrix_E_inv(1,2) = -gamma1 * k / v_alpha
  matrix_E_inv(1,3) = 1.0 / (2. * miu )
  matrix_E_inv(1,4) = k / (2. * miu* v_alpha)
  matrix_E_inv(2,1) = gamma1 * k / v_beta
  matrix_E_inv(2,2) = 1.0
  matrix_E_inv(2,3) = -k / (2. * miu * v_beta)
  matrix_E_inv(2,4) = -1.0 / (2. * miu )
  matrix_E_inv(3,1) = 1.0
  matrix_E_inv(3,2) = -gamma1 * k / v_alpha
  matrix_E_inv(3,3) = -1.0 / (2. * miu)
  matrix_E_inv(3,4) = k / (2. * miu * v_alpha)
  matrix_E_inv(4,1) = -gamma1 * k / v_beta
  matrix_E_inv(4,2) = 1.0
  matrix_E_inv(4,3) = k / (2. * miu* v_beta)
  matrix_E_inv(4,4) = -1.0 / (2. * miu )

  matrix_E_inv = matrix_E_inv * 0.5 * gamma

end subroutine cal_E_inv

subroutine cal_E_inv_par(omega, ray_p,alpha,beta,rho,matrix_E_inv, par_types)
  implicit none
  real(dp), intent(in)                :: ray_P, omega,alpha,beta,rho
  complex(dcp),intent(out)            :: matrix_E_inv(4,4)
  integer(c_int), intent(in)          :: par_types  
  
  ! local
  real(dp)                            ::  k,miu   
  complex(dcp)                        :: v_alpha, v_beta , k_alpha, k_beta,&
                                         gamma, gamma1, gamma3  

  ! elastic modulus
  miu = rho * beta**2 
  k = omega * ray_p
  k_alpha = omega / alpha
  k_beta = omega / beta
  v_alpha = sqrt(k**2 - k_alpha**2)
  v_beta = sqrt(k**2 - k_beta**2)     
  gamma = 2 * k**2 * beta**2 / omega**2
  gamma1 = 1 - 1 / gamma
  gamma3 = 1.0_dp / (gamma - 2_dp)
  
  if (par_types == 3) then
    !beta
    matrix_E_inv(:,:) = (0.0d0,0.0d0)
    matrix_E_inv(1,1) = -1.0
    matrix_E_inv(1,2) = -k / v_alpha
    matrix_E_inv(2,1) = k / v_beta *(1 - gamma1*gamma3)
    matrix_E_inv(2,2) = 1.0
    matrix_E_inv(2,3) = k * gamma3 / 2 / miu / v_beta
    matrix_E_inv(3,1) = 1.0
    matrix_E_inv(3,2) = -k / v_alpha
    matrix_E_inv(4,1) = -k / v_beta *(1 - gamma3*gamma1)
    matrix_E_inv(4,2) = 1.0
    matrix_E_inv(4,3) = -k *gamma3 / 2 / miu / v_beta
    matrix_E_inv = matrix_E_inv * gamma / beta
  elseif(par_types == 1) then
    !rou
    matrix_E_inv(:,:) = (0.0d0,0.0d0)
    matrix_E_inv(1,3) = -1.0
    matrix_E_inv(1,4) = -k / v_alpha
    matrix_E_inv(2,3) = k / v_beta
    matrix_E_inv(2,4) = 1.0
    matrix_E_inv(3,3) = 1.0
    matrix_E_inv(3,4) = -k / v_alpha
    matrix_E_inv(4,3) = -k / v_beta
    matrix_E_inv(4,4) = 1.0
    matrix_E_inv = matrix_E_inv * gamma / 4.0_dp / rho / miu
  elseif(par_types == 2) then
    !alpha
    matrix_E_inv(:,:) = (0.0_dp,0.0_dp)
    matrix_E_inv(1,2) = gamma1
    matrix_E_inv(1,4) = -0.5 / miu
    matrix_E_inv(3,2) = gamma1
    matrix_E_inv(3,4) = -0.5 / miu
    matrix_E_inv = matrix_E_inv / alpha *&
                   k_alpha**2 * k**3 / k_beta**2 / v_alpha**3
  else if(par_types == 4) then 
    ! thick
    matrix_E_inv(:,:) = (0.0_dp,0.0_dp)
  else 

  endif
end

subroutine irfft(inp,out,n)
  !! Inverse FFT 
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
     
  integer(c_int),intent(in)           :: n
  complex(dcp),intent(in)             :: inp(n/2+1)
  real(dp),intent(inout)              :: out(n)

  ! fft local
  type(C_PTR)                         :: plan
  complex(dcp)                        :: in(n/2+1)

  ! save input arrays
  in(:) = inp(:)

  ! choose plan 
  plan = fftw_plan_dft_c2r_1d(n,in,out,FFTW_BACKWARD)

  call fftw_execute_dft_c2r(plan, in, out)
  call fftw_destroy_plan(plan)
  out(:) = out(:) / n
end subroutine irfft

end module RFModule