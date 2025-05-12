module genNoise_mod
!  subroutines needed to generate sythetic time series -- requires a random number generator
!     and some filters (along with convolution)
  implicit none
  
  private
  public :: genNoise, convolv, genWhite, convolv0, convolvT
  
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
  
subroutine DoIt2(n,H,F,O)
use OMP_LIB
use iso_fortran_env
integer, intent(in) :: n
real(kind=real64), intent(in):: H(:),F(:)
real(kind=real64), intent(out) :: O(:)
integer :: i,j,k
!$OMP PARALLEL
!$OMP DO
!    print*,"num of threads",OMP_get_num_threads(),' max num of threads',OMP_get_max_threads(), &
 !   OMP_in_parallel(),OMP_get_thread_num()
do i=1,n
  O(i)=0.0
  do j=1,i
    O(i)=O(i) + H(j)*F(i-j+1)
  end do

!    if (i .eq. n) print*,i,O(i),"num of threads",OMP_get_num_threads(),' max num of threads',OMP_get_max_threads()
!    if (i .eq.n) print*,OMP_in_parallel(),OMP_get_thread_num()

end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine DoIt2
   subroutine convolvT(npts,ilen,H,ran,out)
! time domain convolution 
    use iso_fortran_env
    use OMP_LIB
!    implicit none
    integer :: i,j,k
    real(kind=real64), intent(in) :: ran(:)
    real(kind=real64), intent(out) :: out(:)
    real(kind=real64), intent(in) :: H(:)
    real(kind=real64) :: minH
    integer, intent(in) ::   ilen, npts

! !! split in two --- 
    if (ilen .le. npts ) then
!$OMP PARALLEL DO
      do i=1,npts
        out(i)=0.0
        if (i .lt. ilen) then
          do j=1,i
            out(i)=out(i)+ran(j)*H(i-j+1)
          end do
        else
          do j=1,i
            out(i)=out(i)+ran((i+1-j))*H(j)
          end do
        end if
      end do
!OMP END PARALLEL DO
    else
!$OMP PARALLEL DO
      do i=1,npts
        out(i)=0.0
        do j=1,i
          out(i)=out(i)+ran(j)*H(i-j+1)
        end do
      end do
!OMP END PARALLEL DO
    end if
  end subroutine convolvT
  
  subroutine convolvF(n,ilen,filt,data,out) 
!  frequency domain convolution
!  Attempt to use dfft routines to follow that algorthm in NR program convlv
    use iso_fortran_env
!    implicit none
    integer, intent(in) :: n, ilen
    real(kind=real64), intent(in) :: filt(:),data(:)
    real(kind=real64), intent(out) :: out(:)
    real(kind=real64), allocatable :: filt2(:),data2(:)
    real(kind=real64), allocatable :: wsave(:), freal(:),fimag(:),dreal(:),dimag(:),oreal(:), &
      oimag(:)
    real(kind=real64), allocatable :: out2(:)
    real(kind=real64) :: fzero,dzero,minfilt
    integer :: i,npts,iodd
    real(kind=real64) :: pi
    
    pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062

   
! figure-out the effective length of filter function... when does it go to zero
 


!    print*,' effective length of filter function; needs to be odd',ilen  
    npts=n+ilen
!  force npts to be power of 2
    npts=int(alog(float(npts)) / alog(2.) ) + 1
    npts=int(2**npts)
!    npts=n+ilen
 !   npts=int(npts*1.0)
!    print*,' n ilen npts',n,ilen,npts
    allocate(filt2(npts))
    allocate(data2(npts))
    allocate(out2(npts))
    data2=0.0d+0
    filt2=0.0d+0
    do i=1,n
      data2(i)=data(i)
    end do
!    print*,'ilen n npts sizefilt sizefilt2 sizedata ',ilen,n,npts,size(filt),size(filt2),size(data)
    do i=1,ilen
      filt2(i)=filt(i)
    end do

!  redistribute filt2 stuff per NR program convlv
      do  i=1,(ilen-1)/2
        filt2(npts+1-i)=filt2(ilen+1-i)
      end do
      do i=(ilen +3)/2,npts-(ilen-1)/2
        filt2(i)=0.0
      end do

    allocate(wsave(3*2*npts/2+15))
    call dzffti(npts*2/2,wsave)    !! initialize
!!  fft
    allocate(freal(npts))
    allocate(fimag(npts))
    allocate(dreal(npts))
    allocate(dimag(npts))
    allocate(oreal(npts))
    allocate(oimag(npts))
!    do i=1,npts
!      write(36,*)i,data2(i),filt2(i)
!    end do
 
    call dzfftf(2*npts/2,filt2,fzero,freal,fimag,wsave)       
    call dzfftf(2*npts/2,data2,dzero,dreal,dimag,wsave)   
    print*,'dzero and fzero', dzero,fzero   
    oreal=dreal*freal-dimag*fimag
    oimag=freal*dimag+fimag*dreal
!  inverse fft
!    allocate(out2(2*npts))
!    fzero=1.0d+0/((float(n)/2.5)**2)
!    fzero=fzero*dzero
    fzero=0.0d+0
!    fzero=dzero
    call dzfftb(2*npts/2,out2,fzero,oreal,oimag,wsave)
    do i=1,n
!        out(i)=sqrt(2.0)*out2(i)*float(n)
!        out(i)=out2(i)
!        out(i)=out2(i)*(dble(float(n)))/(dsqrt(pi))
        out(i)=out2(i)*(dble(float(npts/2)))      !-dzero
!        out(i)=out2(i)
    end do
    deallocate(wsave)
    deallocate(filt2)
    deallocate(data2)
    deallocate(freal)
    deallocate(fimag)
    deallocate(dreal)
    deallocate(dimag)
    deallocate(oreal) 
    deallocate(oimag)
  end subroutine 
  subroutine convolvTx(npts,H,ran,ilen)
! time domain convolution 
    use iso_fortran_env
!    implicit none
    integer :: npts,i,j,ilen,k
    real(kind=real64), intent(inout) :: ran(:)
    real(kind=real64) :: ranin(size(ran))
    real(kind=real64), intent(in) :: H(:)
    real(kind=real64) :: minH
! figure-out the effective length of filter function... when does it go to zero
!    print*,'ilen=',ilen
    ranin=ran
    do concurrent (i=1:npts)
!    do i=1,npts
      ran(i)=0.0
!      do  j=1,i
 !     k=i
!      if (i .gt. ilen) k=ilen
      if (i .lt. ilen) then
        do concurrent (j=1:i)
          ran(i)=ran(i)+ranin(j)*H(i-j+1)
        end do
      else
        do concurrent (j=1:ilen)
          ran(i)=ran(i)+ranin((i+1-j))*H(j)
        end do  
      end if    
    end do
  end subroutine convolvTx
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
  
  subroutine convolv1(n,filt,data,out)
!  frequency domain convolution
    use iso_fortran_env
!    implicit none
    integer, intent(in) :: n
    real(kind=real64), intent(in) :: filt(:),data(:)
    real(kind=real64), intent(out) :: out(:)
    real(kind=real64), allocatable :: filt2(:),data2(:)
!    real(kind=real64), intent(out) ;; out(:)
    real(kind=real64), allocatable :: wsave(:), freal(:),fimag(:),dreal(:),dimag(:),oreal(:), &
      oimag(:),out2(:)
    real(kind=real64) :: fzero,dzero
    integer :: i,npts,ndif
    real(kind=real32) :: pi
    pi=3.14159265
!  force npts be a power of 2  !! doesn't have to be....
      npts=n
      npts=int(alog(float(npts))/alog(2.))+1  
 
      npts=1*(2**npts)  
      
      allocate(filt2(2*npts))
      allocate(data2(2*npts))

!!  zero pad the inputs to the fft 
      data2=0.0d+0
      filt2=0.0d+0
      ndif=npts-n
!      print*,' in convolv1 npts, n, ndif  size data2 size data',npts,n,ndif,size(data2),size(data)
!      print*,1+npts/2-1+ndif/2,n+npts/2-1+ndif/2
      do i=1,n
        data2(i+npts/2-1+ndif/2)=data(i)
        filt2(i+npts/2-1+ndif/2)=filt(i)
!        if (i .ge. 9*n/10) then     !!!  cosine taper last 10% of filter
!           filt2(i+npts/2-1+ndif/2)=filt2(i+npts/2-1+ndif/2)*cos((pi/2.0)*(i-9*n/10)/(n/10))
!           data2(i+npts/2-1+ndif/2)=data2(i+npts/2-1+ndif/2)*cos((pi/2.0)*(i-9*n/10)/(n/10))
!        end if
      end do

      allocate(wsave(3*2*npts+15))
      call dzffti(npts*2,wsave)    !! initialize
!!  fft
      allocate(freal(npts))
      allocate(fimag(npts))
      allocate(dreal(npts))
      allocate(dimag(npts))
      allocate(oreal(npts))
      allocate(oimag(npts))
 
      call dzfftf(2*npts,filt2,fzero,freal,fimag,wsave)       
      call dzfftf(2*npts,data2,dzero,dreal,dimag,wsave) 
      oreal=dreal*freal-dimag*fimag
      oimag=freal*dimag+fimag*dreal
!      do i=1,npts
!          j=i
!        write(78,*)i,j,(oreal(j)**2+oimag(j)**2), dlog(oreal(j)**2+oimag(j)**2)/alog(10.),  &
!        dlog(freal(j)**2+fimag(j)**2)/alog(10.),dlog(dreal(j)**2+dimag(j)**2)/alog(10.)
!        write(76,*)i,freal(i),fimag(i),dlog(freal(j)**2+fimag(j)**2)/alog(10.)
!      end do
!  inverse fft
      allocate(out2(2*npts))

      dzero=1.0d+0/((float(n)/2.5)**2)     
      call dzfftb(2*npts,out2,dzero,oreal,oimag,wsave) 
!      do i=1,2*npts
!        write(70,*)i,out2(i)
!      end do
      do i=1,n
        out(i)=dsqrt(2.0d+0)*out2(npts+i+ndif/1-2)*float(n)
      end do
      deallocate(wsave)
      deallocate(filt2)
      deallocate(data2)
      deallocate(freal)
      deallocate(fimag)
      deallocate(dreal)
      deallocate(dimag)
      deallocate(oreal)
      deallocate(oimag)
      
  end subroutine convolv1
 
  subroutine convolv4(n,filt,data,out)
!  frequency domain convolution
!  Attempt to use dfft routines to follow that algorthm in NR program convlv
    use iso_fortran_env
!    implicit none
    integer, intent(in) :: n
    real(kind=real64), intent(in) :: filt(:),data(:)
    real(kind=real64), intent(out) :: out(:)
    real(kind=real64), allocatable :: filt2(:),data2(:)
    real(kind=real64), allocatable :: wsave(:), freal(:),fimag(:),dreal(:),dimag(:),oreal(:), &
      oimag(:)
    real(kind=real64), allocatable :: out2(:)
    real(kind=real64) :: fzero,dzero,minfilt
    integer :: i,npts,ilen,iodd
    real(kind=real64) :: pi
    
    pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062

   
! figure-out the effective length of filter function... when does it go to zero
    minfilt=maxval(dsqrt(filt**2))*1.0e-04
    ilen=n-100
    do i=1,n-100
     
       if (dsqrt(filt(n-100+1-i)**2) .ge. minfilt) then
         ilen=n-100+1-i
         exit
       end if
    end do
    ilen=ilen*2
    if (ilen .gt. n)ilen=n-100
    iodd=int(ilen/2)
    iodd=ilen-2*iodd
    if (iodd .eq. 0) ilen=ilen+1
    if (ilen .gt. n) ilen=ilen-2

!    print*,' effective length of filter function; needs to be odd',ilen  
    npts=n+ilen
!  force npts to be power of 2
    npts=int(alog(float(npts)) / alog(2.) ) + 1
    npts=int(2**npts)
!    npts=n+ilen
 !   npts=int(npts*1.0)
!    print*,' n ilen npts',n,ilen,npts
    allocate(filt2(npts))
    allocate(data2(npts))
    allocate(out2(npts))
    data2=0.0d+0
    filt2=0.0d+0
    do i=1,n
      data2(i)=data(i)
    end do
!    print*,'ilen n npts sizefilt sizefilt2',ilen,n,npts,size(filt),size(filt2)
    do i=1,ilen
      filt2(i)=filt(i)
    end do

!  redistribute filt2 stuff per NR program convlv
      do  i=1,(ilen-1)/2
        filt2(npts+1-i)=filt2(ilen+1-i)
      end do
      do i=(ilen +3)/2,npts-(ilen-1)/2
        filt2(i)=0.0
      end do

    allocate(wsave(3*2*npts/2+15))
    call dzffti(npts*2/2,wsave)    !! initialize
!!  fft
    allocate(freal(npts))
    allocate(fimag(npts))
    allocate(dreal(npts))
    allocate(dimag(npts))
    allocate(oreal(npts))
    allocate(oimag(npts))
    do i=1,npts
      write(36,*)i,data2(i),filt2(i)
    end do
 
    call dzfftf(2*npts/2,filt2,fzero,freal,fimag,wsave)       
    call dzfftf(2*npts/2,data2,dzero,dreal,dimag,wsave)   
    
    oreal=dreal*freal-dimag*fimag
    oimag=freal*dimag+fimag*dreal
!  inverse fft
!    allocate(out2(2*npts))
!    fzero=1.0d+0/((float(n)/2.5)**2)
!    fzero=fzero*dzero
    fzero=0.0d+0
    call dzfftb(2*npts/2,out2,fzero,oreal,oimag,wsave)
    do i=1,n
!        out(i)=sqrt(2.0)*out2(i)*float(n)
!        out(i)=out2(i)
!        out(i)=out2(i)*(dble(float(n)))/(dsqrt(pi))
        out(i)=out2(i)*(dble(float(npts/2)))+dzero
!        out(i)=out2(i)
    end do
    deallocate(wsave)
    deallocate(filt2)
    deallocate(data2)
    deallocate(freal)
    deallocate(fimag)
    deallocate(dreal)
    deallocate(dimag)
    deallocate(oreal) 
    deallocate(oimag)
  end subroutine convolv4

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
