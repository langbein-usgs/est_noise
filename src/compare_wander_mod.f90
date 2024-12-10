module compare_wander_mod

  use iso_fortran_env
  use genNoise_mod, only: convolv0
  implicit none
  private
  real(kind=real32), public :: exp1,amp1,alpha,exp2,amp2,sig,fl,fh,amp_bp,fmiss
  integer, public :: np
  real(kind=real64), allocatable, public :: wh(:),pl1(:),pl2(:),bp(:), t(:)  
  public :: wander, MakeNoise, chron
  contains
  
  subroutine wander(tim_wan,wan_data,ntt,dat,ix,dt_sam)

!  calculates the wander for specificed intervals
!   wan_data(T)=sum (dat(T+t)-dat(t))^2 / N
!  tim_wan are the time interval (days) to estimate wander--input
!  wan_dat are the estimates of wander for specified tim_wan--output
!  dat observations (missing data 99999)
!  max_wan and max_flt  array dimensions
!  ntt  number of periods in tim_wan
!  ic number of data
!  dt_sam sampling interval in days

  real(kind=real32), intent(in) :: tim_wan(:),dat(:)
  real(kind=real32), intent(out) :: wan_data(:)
  real(kind=real64), intent(in) :: dt_sam
  integer, intent(in) :: ntt
  real(kind=real32) :: fmiss,tau
  integer :: j,itau,itau_min,itau_max,k_tau,n,ix,i 
  integer :: max_flt,max_wan
!  missing data is 99999.
  fmiss=99999.
  do  j=1,ntt
    wan_data(j)=0.0
    n=0
    tau=tim_wan(j)
    itau=int(tau/dt_sam)
    if (j .eq. 1) itau_min=int(tim_wan(1)/dt_sam)
    if (j .gt. 1) itau_min=int((tim_wan(j-1))/dt_sam +1)
    itau_max=int(tim_wan(j)/dt_sam)
!    print*,tim_wan(j),itau_min,itau_max
    do  k_tau=itau_min,itau_max
      do  i=1,ix-k_tau
        if ((dat(i) .ne. fmiss ) .and. (dat(i+k_tau) .ne. fmiss)) then
           n=n+1
           wan_data(j)=wan_data(j)+(dat(i+k_tau)-dat(i))**2
         end if
      end do
    end do
    if (n .ne. 0) wan_data(j)=(wan_data(j)/float(n))
    if (n .eq. 0) wan_data(j)=fmiss
!     print*,j,tim_wan(j),itau,itau_min,itau_max,n,wan_data(j)

  end do
  end subroutine
  subroutine MakeNoise(ran,Modtype,npts,dt_sam,iseed)
 
    real(kind=real32), intent(out) :: ran(:)
    real(kind=real64), intent(in) :: dt_sam
    real(kind=real64), allocatable :: ran1(:), ran2(:), ran3(:), ran4(:),filt(:)
    integer, intent(in) :: npts,iseed
    character(len=1) :: Modtype
    real(kind=real64) :: dts, ZBQLNOR
    real(kind=real32) :: ts,sum
    integer :: i
    allocate(ran1(npts+100))
    allocate(ran2(npts+100))
    allocate(ran3(npts+100))
    allocate(ran4(npts+100))
    allocate(filt(npts+100))
    ts=sngl(dt_sam/365.25)
    dts=(dt_sam/365.25)
    ran1=0.0d+0
    ran2=0.0d+0
    ran3=0.0d+0
    ran4=0.0d+0
    filt=0.0d+0
    9     continue
    call srand(iseed)
    if ((ModType .eq. 'n') .or. (ModType .eq. 'q')) then
!   quadrature noise
      if (sig .ne. 0.0) then
        do i=1,npts+100
          ran1(i)=sig*sngl(ZBQLNOR(0.0d+0,1.0d+0))
        end do
        call convolv0(npts+100,wh,ran1)
      end if
      if (amp1 .ne. 0.0) then
        do  i=1,npts+100
          ran2(i)=amp1*sngl(ZBQLNOR(0.0d+0,1.0d+0))*((ts)**(exp1/4))
        end do
        call convolv0(npts+100,pl1,ran2)
      end if
      if (amp2 .ne. 0.0) then
        do  i=1,npts+100
          ran3(i)=amp2*sngl(ZBQLNOR(0.0d+0,1.0d+0))*((ts)**(exp2/4))
        end do
        call convolv0(npts+100,pl2,ran3)
      end if

      if (amp_bp .ne. 0.0) then
        do  i=1,npts+100
          ran4(i)=amp_bp*sngl(ZBQLNOR(0.0d+0,1.0d+0))
        end do
        call convolv0(npts+100,bp,ran4)
      end if
      filt=ran1+ran2+ran3+ran4

      ran1=filt 

    else
! additive noise
      do i=1,npts+100
        ran1(i)=1.0*(ZBQLNOR(0.0d+0,1.0d+0))
        filt(i)=wh(i)*sig+amp1*pl1(i)*((ts)**(exp1/4))+amp2*pl2(i)*((ts)**(exp2/4))+amp_bp*bp(i)        
      end do
      call convolv0(npts+100,filt,ran1)
    end if

!   Resample the noise at time of measurement
   do i=1,npts
     if (ran(i) .ne. fmiss) ran(i)=ran1(i+100)

   end do

     
  

    deallocate(filt)
    deallocate(ran1)
    deallocate(ran2)
    deallocate(ran3)
    deallocate(ran4)
  end subroutine
  
  subroutine chron(t,ic)
!  puts data into chronological order
!  code riped-off from numerical recipes piksrt.f
     real(kind=real32) :: t(:),tx
     integer :: i,j,ic
      do  j=2,ic
         tx=t(j)

         do  i=j-1,1,-1
           if (t(i) .eq. tx) t(i)=t(i)+0.001
           if (t(i) .lt. tx) go to 10
           t(i+1)=t(i)
         end do
         i=0
10       t(i+1)=tx
      end do
      return
      end

end module 


