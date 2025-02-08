module genNoise_mod
!  subroutines needed to generate sythetic time series -- requires a random number generator
!     and some filters (along with convolution)
  implicit none
  
  private
  public :: genNoise, convolv, genWhite, convolv0
  
  contains 
!
  subroutine genNoise(ModType,t,d,iseed,dt)

!  substitutes noise for real data
!
!  ModType, n, a, c, f
!        n and q is standard, noise summed from independent random number generators
!         a, c, or f sums the transfer function and uses single random number generator
!  t is time
!  d is data-  which will be filled with noise
!  ic and max_data are number of data and the dimensions of arrays 
  use iso_fortran_env
  use filterfunc_mod
  character(len=1) ,intent(in):: ModType
  integer, intent(inout) :: iseed
  real(kind=real64), intent(in) :: t(:), dt
  real, intent(in out) :: d(:)
  real :: expon1,expon2,amp1,amp2,alpha,sig,fl,fh,f_amp,f_tem,ts
  real, allocatable ::  ran(:,:), ran0(:)
  real(kind=real64), allocatable :: filt(:,:),temp(:),temp2(:)
  complex(kind=real64), allocatable :: ran1(:)
!  real(kind=real64), allocatable :: f1(:), f2(:), f3(:), f4(:)
  integer :: np, nmax,npts, i,j
  integer, parameter :: ispinUp=100    !  number of ghost data before the time-series starts
!  time declarations
  real :: time0,dtime,ddtime0
!!  input the parameters for the noise model to be created
  print*,' Program creates colored noise with'
  print*,'  PSD   Po / ((alpha^n + f ^ n)) Plus white noise'
  print*,'     plus more power-law noise of P1/f^j'
  print*,'     pluse band-passed filterd noise'
  print*,' '
  print*,' First power law/Gauss-Markov parameters'
  print*,' input the exponent n'
  read(5,*)expon1
  write(25,fmt="(f15.3,'     # input noise 1 exp term')")expon1
  print*,' input term proportional to Po'
  read(5,*)amp1
  write(25,fmt="(f15.3,'     # input noise 1 amp. term')")amp1
  print*,' input the time constant alpha in radians/year'
  read(5,*)alpha
  write(25,fmt="(f15.5,'     # input noise 1 GM term')")alpha
  print*,' '
  print*,' Second power law parameters'
  print*,' input the exponent n'
  read(5,*)expon2
  write(25,fmt="(f15.3,'     # input noise 2 exp term')")expon2
  print*,' input term proportional to Po'
  read(5,*)amp2
  write(25,fmt="(f15.3,'     # input noise 2 amp. term')")amp2
  print*,' '
  print*,' Input the white noise'
  read(5,*)sig
  write(25,fmt="(f15.3,'     # input white noise')")sig
  print*,' '
  print*,' Creates a second set of noise consisting of band-passed'
  print*,'  band passed filtered, white noise'
  print*,' Input the lower and upper limits of the filter in c/yr'
  read(5,*)fl,fh
  if (fl .gt. fh) then
     f_tem=fl
     fl=fh
     fh=f_tem
  end if
  write(25,fmt="(2f8.2, '  #lower and upper limits of the BP filter')") fl,fh
  print*,' Input the number of pole of the filter, 1,2,..'
  read(5,fmt="(i3)") np
  write(25,fmt="(i3, '  #  Num of poles')")np
  print*,' Input the amplitude of the noise'
  read(5,*)f_amp
  write(25,fmt="(f15.3,'     # Input Band Pass amp')") f_amp 
  


  ts=(dt/365.25)
  npts=size(t)
  nmax=size(t)+ispinUp

  allocate(ran0(nmax))
  allocate(ran(4,nmax))
  allocate(filt(4,nmax))
  allocate(temp(nmax))
  allocate(temp2(nmax))

  ran=0.0
  filt=0.0
  ran0=0.
  temp=0.0
  temp2=0.0

  if (sig .ne. 0) then  
!! white noise impulse response
    filt(1,1)=1.0d+0
    do i=2,nmax
      filt(1,i)=0.0d+0
    end do
   end if
   
  
  if (amp1 .ne. 0) then
    temp=0.0
!  First power law impulse response
     call frac_diff(temp,expon1,alpha,ts,nmax)
     filt(2,:)=temp
   end if
  if (amp2 .ne. 0) then
!  second power law impulse response
     temp=0.0
     call frac_diff(temp,expon2,0.0,ts,nmax)
     filt(3,:)=temp
   end if
   if (f_amp .ne. 0.0) then
!! bandpassed filter
     temp=0.0
     call band_pass_filt(dble(ts),fl,fh,np,nmax,temp,1.0)
     filt(4,:)=temp
   end if

!
!  Depending upon noise model type, either add all of the implus responses
!    and convolve with noise, or convolve noise with each impulse response, and then add
!
  print*,' ModType ',ModType

  if ((ModType .eq. 'n' ) .or. (ModType .eq. 'q' )) then
    print*,' standard noise of 4 independent noise models'
    if (sig .ne. 0) then
      call genWhite(nmax,iseed,ran0)
      ran(1,:)=sig*ran0
      print*, "sig=",sig
      temp=filt(1,:)
      temp2=dble(ran(1,:))
      call convolv0(nmax,temp,temp2)
      ran(1,:)=sngl(temp2)
!      call convolv(nmax,filt(1,:),ran(1,:))
    end if
    if (amp1 .ne. 0) then
       call genWhite(nmax,iseed,ran0)
       ran(2,:)=amp1*ran0*((ts)**(expon1/4))
       print*,' amp1 ggm',amp1,alpha
       temp=filt(2,:)
       temp2=dble(ran(2,:))
       call convolv0(nmax,temp,temp2)
       ran(2,:)=sngl(temp2)
!       call convolv(nmax,filt(2,:),ran(2,:))
    end if
    if (amp2 .ne. 0) then
       call genWhite(nmax,iseed,ran0)   
       ran(3,:)=amp2*ran0*((ts)**(expon2/4))
       print*, 'amp2=',amp2,((ts)**(expon2/4))
       temp=filt(3,:)
       temp2=dble(ran(3,:))
       call convolv0(nmax,temp,temp2)
       ran(3,:)=sngl(temp2)
!       call convolv(nmax,filt(3,:),ran(3,:))
    end if
    if (f_amp .ne. 0) then 
       call genWhite(nmax,iseed,ran0) 
       print*,'bp_amp' ,f_amp 
       ran(4,:)=f_amp*ran0
       temp=filt(4,:)
       temp2=dble(ran(4,:))
       call convolv0(nmax,temp,temp2)
       ran(4,:)=sngl(temp2)
       call convolv(nmax,filt(4,:),ran(4,:))
    end if
    
!    print*,'start sumation'
    do  i=1,nmax
      ran(1,i)=ran(1,i)  +ran(2,i)+ran(3,i)+ran(4,i)
    end do
!    print*,' end sumation'
  else
  
     print*,' sum impulse responses and then add'
     do i=1,nmax
        filt(1,i)=filt(1,i)*sig+amp1*filt(2,i)*((ts)**(expon1/4))+ &
           amp2*filt(3,i)*((ts)**(expon2/4))+f_amp*filt(4,i)
     end do

      ran0=ran(1,:)
      call genWhite(nmax,iseed,ran0) 
      ran(1,:)=ran0

     
!     call convolv(nmax,filt(1,:),ran(1,:))
     temp=filt(1,:)
     temp2=dble(ran(1,:))
     call convolv0(nmax,temp,temp2)
     ran(1,:)=sngl(temp2)
  end if
  
  
! grand total
  do  i=1,npts
     d(i)=ran(1,i+ispinUp)
!     write(16,*)i,ran(1,1+ispinUp),real(ran1(i+ispinUp)),aimag(ran1(i+ispinUp))
  end do

  
  deallocate(ran0)
  deallocate(ran)
  deallocate(filt)
  deallocate(temp)
  deallocate(temp2)
  end subroutine genNoise 
  
  subroutine convolv0(npts,H,ran)
! time domain convolution 
    use iso_fortran_env
    implicit none
    integer :: npts,i,j
    real(kind=real64), intent(inout) :: ran(:)
    real(kind=real64) :: ranin(size(ran))
    real(kind=real64), intent(in) :: H(:)

    ranin=ran
    do concurrent (i=1:npts)
!    do i=1,npts
      ran(i)=0.0
!      do  j=1,i
      do concurrent (j=1:i)
        ran(i)=ran(i)+ranin(j)*H(i-j+1)
      end do
    end do
  end subroutine convolv0
  

  
  subroutine convolv(npts,H,ran)
!  fft version of convolution... faster than time domain convolution
    use iso_fortran_env
    implicit none
    integer :: npts,i,j,kmax,len
    real, intent(inout) :: ran(:)
    real(kind=real64), intent(in) :: H(:)

    real(kind=real64), allocatable :: wsave(:)
    real(kind=real64), allocatable :: a1(:),b1(:), a2(:),b2(:),ranin(:),hh(:)
    real(kind=real64) :: azero
    complex(kind=real64) :: hc,rc,hr


 
!  make length of time series a power of two and pad excess with zero
    len=npts
    len=int((log(float(len))/log(2.))+0)+2
    len=2**len
    print*,'npts=',npts,' len=',len
    kmax=(len)/2+1
    allocate(wsave(3*len+15))
    allocate(a1(kmax))
    allocate(b1(kmax))
    allocate(a2(kmax))
    allocate(b2(kmax))
    allocate(hh(len))
    allocate(ranin(len))

!  zero fill ends  (ie, pad time series)
!
    do i=1,len
      hh(i)=0.00
      ranin(i)=0.0
    end do
!    do i=1,(len-npts)/2
!       hh(i)=0.0
!       if (npts+1-i .gt. 0 ) hh(npts+1-i)=0.0
!       ranin(i)=0.0
!       if (npts+1-i .gt. 0 ) ranin(npts+1-i)=0.0
!    end do
    do i=1,npts
      hh(i+(len-npts)/2)=H(i)
      ranin(i+(len-npts)/2)=dble(ran(i))
    end do
!  fft of both random noise and the filter/impulse response
    call dzffti(len,wsave)
    call dzfftf(len,hh,azero,a2,b2,wsave)
!    print*,'fft of filter',azero 
    call dzffti(len,wsave)
    call dzfftf(len,ranin,azero,a1,b1,wsave)
!    print*,'fft of white noise',azero

!  multiply the two fft 
    do i=1,kmax
!       write(14,*)i,sqrt(a1(i)**2+b1(i)**2),(180./3.14)*atan2(b1(i),a2(i)), &
!                    sqrt(a2(i)**2+b2(i)**2),(180./3.14)*atan2(b2(i),a2(i))
       rc=cmplx(a1(i),b1(i))
       hc=cmplx(a2(i),b2(i))
       hr=hc*rc   !  the multiplication  of two complex numbers for convolution
!       a1(i)=(a1(i)*a2(i) - b1(i)*b2(i))*float(len/2)   ! alternative multiplication
!       b1(i)=(b1(i)*a2(i) + b2(i)*a1(i))*float(len/2)


       a1(i)=real(hr)*float(len/2)
       b1(i)=(-1.0d+0)*aimag(hr)*float(len/2)
    end do
!  inverse fft
    call dzffti(len,wsave)
    call dzfftb(len,ranin,0.0d+0,a1,b1,wsave) 
!                   (n,r,azero,a,b,wsave
!    print*,' invfft ',azero

!  output the convolution accounting for the padded time series
    do i=1,npts
      ran(i)=sngl(ranin(i+(len-npts)/2)) 
!      write(16,*) i-100,ranin(i),ran(i)
    end do
!    ran=sngl(ranin)*sqrt(1.0)
    deallocate(wsave)
    deallocate(a1)
    deallocate(b1)
    deallocate(a2)
    deallocate(b2)
    deallocate(hh)
    deallocate(ranin)
  end subroutine convolv
  
  subroutine genWhite(npts,iseed,ran)
!  generates white noise with ave=0.0 and stan dev = 1
!  serves as a wrapper for ZBQLNOR as that routing has a bug
!   ZBQLNOR will occasionally output a Nan
    use iso_fortran_env
    implicit none
    integer :: i  
    integer, intent(in) :: npts,iseed
    real(kind=real32), intent(out)  :: ran(npts)
    real(kind=real64) ::  ZBQLNOR
    call srand(iseed)  
    do i=1,npts
      ran(i)=ZBQLNOR(0.0d+0,1.0d+0)
!  is the noise within bounds; if not repeat
      if ((ran(i) > 1000.d+0) .or. (ran(i) < -1000.d+0)) then
        ran(i)=sngl(ZBQLNOR(0.0d+0,1.0d+0))
      end if
!  is the noise within bounds; if not repeat
      if ((ran(i) > 1000.d+0) .or. (ran(i) < -1000.d+0)) then
        ran(i)=sngl(ZBQLNOR(0.0d+0,1.0d+0))
      end if
    end do
  end subroutine genWhite
  
end module genNoise_mod
