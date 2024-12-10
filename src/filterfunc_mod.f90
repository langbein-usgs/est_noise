module filterfunc_mod
!  filters use to create synthetic noise and/or to create data covariance matrices
  use iso_fortran_env
  use alloc_mod
  implicit none
  private
 
  public :: frac_diff, band_pass_filt, pow_law_cov,band_pass_cov, makeMatrix
  
  contains 
  
  subroutine frac_diff(crr,expon,alpha,ts,npts)
! use Hosking routine to create fractional differences
!  crr are the fractional difference
!  expon is the exponent
!  alpha is the GM freq in radians/yr
!  ts is the sampling interval in years
!  npts is length of time series
!  nmax is dimension of crr
!  Hosking, J. R. M. (1981). Fractional Differencing. Biometrika, 68(1), 
!     165–176. https://doi.org/10.2307/2335817
!!
    integer, intent(in):: npts
    real, intent(in) :: expon,alpha, ts
    real(kind=real64), intent(out) :: crr(npts)
    real(kind=real64), parameter :: small= 1.0d-10
    real :: ddexp,fracti,fracti_last
    integer :: i,j
    ddexp=0.5*expon
      do  j=1,npts
         if (j .eq. 1) fracti=1.0
         if (j .eq. 2) then
            fracti=ddexp
            fracti_last=fracti
         end if
         if (j .ge. 3) then
            fracti=fracti_last*(ddexp+float(j-2))/float(j-1)
            fracti_last=fracti
         end if
         crr(j)=fracti*dexp(-dble(alpha)*ts*float(j-1))
         if (crr(j) .lt. small) crr(j)=small
      end do
  end subroutine frac_diff  
  
  subroutine band_pass_filt(ts,fl,fh,np,nmax,Fi,amp)
! create a bandpass filter in the time domain
!  define filter in frequency domain then take inverse fft of filter
!    tranfer function to get time domain impulse responce
!  ts is sampling interval in years
!  fl and fh are low and high freq break pts in c/yr
!  np is number of poles 1 to 4
!  nmax is the dimension
!  npts is the number of points in time series
!  F on output is the impulse responce
!  amp is the amplitude of filter.
!
!  Modified May 19, 2015 to use fftpack from netlib.org instead of
!   realft from Numerical Recipes
!
!  Modify Aug 2020 to chose whether to use the invFFT or convolution
!   to generate the BP filter function.
!  Convolution seems to smooth-out the jitter seen in the BP impulse response
!   from the invFFT seen in the time-domain
!
!  for now, the choice is hard wired; (if convolution is deemed better,
!    then it should be possible to compile without the FFT (and invFFT) routines)
!
!  November 2004;
!   changed normalization of band_pass_filt such that sums of squares of filter coefficient  equals 1.00
!   changed band_pass_filter from verson 7 such that I needed to divide fl and fh by 2

    real(kind=real64), intent(in) :: ts
    real, intent(inout) :: fl,fh !,amp
    integer, intent(in) :: np,nmax
    real(kind=real64), intent(out) :: Fi(nmax)
    real, intent(in) :: amp
    complex(kind=real64) :: h1,h2,h3,h4
    complex(kind=real64), allocatable :: hc(:),Fii(:)
    real(kind=real64), allocatable :: a(:),b(:) ,wsave(:)       !,c,d,ddd
    real(kind=real64) :: sumF
    real :: fn, peak,f,f_peak
    integer :: len,i,nf
    integer :: nptt,npoint, nsum,j,k,imax, imin  ! parameters for a smoother
    real :: fmax,fmin,rng1

    print*,'In BP_filt v3',ts,fl,fh,np,nmax,amp
    
    len=nmax
    fh=fh*ts/2.0
    fl=fl*ts/2.0
!    fl=fl/2.0
!    fh=fh/2.0
 !   amp=amp*2.0
    fn=0.5
    nf=int(len/2)
    nf=int((log(float(len/2))/log(2.))+0)
    nf=nf+1
    nf=2**nf
    print*,' Number of frequencies ',nf
    len=2*nf
    print*,' length of time series for band-pass noise is ',len 
    
!!  Allocate arrays
    allocate(hc(len))
    allocate(a(len))
    allocate(Fii(len))
    allocate(b(len))
    allocate(wsave(3*len+15))
!  create the filter function of a bandpass filter in frequency domain
 
    hc(1)=cmplx(0.0d+0,0.0d+0)
    peak=-99999.
    do i=1,nf
      f=(1/float(len))*(i)
       h1=cmplx(1.0,0.0)
       h2=cmplx(0.,f/fl)
       h3=cmplx(1.0,f/fl)
       h4=cmplx(1.,f/fh)
       hc(i)=(h1*h2/(h3*h4))**np

       a(i)=real(hc(i))
       b(i)=-1.0*aimag(hc(i))
       if (sqrt(a(i)**2+b(i)**2) .gt. peak) then
          peak=sqrt(a(i)**2+b(i)**2)
          f_peak=f/ts
       end if
    end do  

      
!  Initialize
    call dzffti(len,wsave)  
!  dfft
    call dzfftb(len,Fii,0.0d+0,a,b,wsave)
!    do  i=1,len
    do i=1,nmax
      Fi(i)=((amp/peak)*Fii(i)*2.0/float(len))*sqrt(0.5/ts)*2.0
    end do
    
    deallocate(hc)
    deallocate(a)
    deallocate(b)
    deallocate(Fii)
 
      
!
!  smooth Fi output using npoint smoother
!
    npoint=3
    nptt=int(npoint/2)
    do i=1,nmax
        wsave(i)=Fi(i)
    end do
    do i=1,nmax
      sumF=0.0d+0
      nsum=0
     do j=1,npoint
        k=i-nptt+j-1
        if (( k .ge. 1) .and. (k .le. nmax)) then
          sumF=sumF+wsave(k)
          nsum=nsum+1
         end if
       end do
 
       if (nsum > 0) Fi(i)=sumF/float(nsum)
 
    end do
   fmax=-99999.
    fmin=999999. 
    do i=1,nmax
        if (Fi(i) .gt. fmax) then
          fmax=Fi(i)
          imax=i
        end if
        if (Fi(i) .lt. fmin) then
          fmin=Fi(i)
          imin=i
        end if
    end do
    rng1=fmax-fmin
    print*,' invFFT range, max, min, point;',rng1,fmax,imax,fmin,imin 
    sumF=0.0
    do i=1,nmax
       sumF=sumF+Fi(i)**2  
    end do
    print*,' Variance of BP filter noise', sumF,' normalize by ',1/sqrt(sumF)
!  normalize to variance such that output variance is 1.0
    Fi=Fi/sqrt(sumF) 
 
         
    deallocate(wsave) 
    fh=2.0*fh/ts
    fl=2.0*fl/ts 

        
  end subroutine band_pass_filt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pow_law_cov(corran,jmax,sig1,sig2,sig3,sig4,covar,mmax,md,t_small)
!  Power Law covariance
    integer, intent(in) :: jmax,mmax,md
    real(kind=real64), intent(out) :: covar(md,md)
    real(kind=real64), intent(in)  :: t_small, corran(md)
    real(kind=real64) :: dcr(mmax), dsum,sumlast,dmax,dmin,dmin0,dcmax,dcmin,dsmall,dtemp
    real(kind=real32), intent(in) :: sig1,sig2,sig3
    real(kind=real32), intent(inout) :: sig4
    real(kind=real32) :: sig4min, const
 
    integer :: irow(md),kmax,kkmax,i,imax,j,ix
!    print*,'in pow_law_cov'

    do  i=1,jmax
    do  j=1,i
        covar(i,j)=0.d+0
        covar(j,i)=0.d+0
    end do
    end do




!   power law covariance
! assumes corran(1) .eq. 0
!  colored noise for exponent
!
! USE FRACTIONAL DIFFERENCE
!       print*,'  jmax=', jmax,' corran(jmax)= ',corran(jmax)-corran(1)
    kmax=int((corran(jmax)-corran(1))/t_small + 0.5)+1
    if (kmax .gt. mmax) then
         print*,' kmax = ',kmax,' exceeds mmax of ',mmax
         stop
    end if
!       time0=secnds(0.0)
!  modify to fix GGM problem (Apr 2021)
    kkmax=kmax
    irowOffset=0
    if (sig4 .ne. 0.0) then
!  put limit on sig4 (GM) freq
      sig4min=0.1*(6.28/corran(jmax))
      if (sig4 .lt. sig4min) sig4=sig4min
       kkmax=kmax+ifix(3*1/(sngl(t_small)*sig4))
       if (kkmax .ge. mmax) kkmax=mmax
      irowOffset=kkmax-kmax
    end if   !!! end GGM fix

    call  frac_diff(dcr,sig3,sig4,sngl(t_small),kkmax)
                
    do  i=mmax,2,-1
         dcr(i)=dcr(i-1)
    end do
    dcr(1)=0

  
    const=(sig2**2)*(t_small**(sig3/2.0))


    do i=1,jmax
        ix=int((corran(i)-corran(1))/t_small+0.5)+1

        irow(i)=ix+irowOffset
    end do
    imax=int((corran(jmax)-corran(1))/t_small+0.5)+1+irowOffset



    call makeMatrix(irow,imax,(dcr),covar,jmax,0,mmax,md)
            

    do  i=1,jmax
      do  j=1,i
        covar(i,j)=covar(i,j)*const
        if (i .eq. j) covar(i,j)=covar(i,j)+sig1**2
        covar(j,i)=covar(i,j)
       end do
    end do
   
 end subroutine pow_law_cov
      
 subroutine band_pass_cov(corran,jmax,t_small,flow,fhigh,npole,mmax,md,sig5,covar)
    integer :: npole, jmax, mmax,md
    real(kind=real64), intent(in) :: t_small, corran(md)
    real(kind=real64), intent(out) :: covar(md,md)
    real(kind=real32) :: flow,fhigh,sig5
    real(kind=real64) :: dsum,sumlast,cr(mmax),time0,timex
    integer :: len, i,j,imax,kk,kmax,k
!!       double precision corran(md),covar(md,md),cr(mmax),t_small
    len=int((corran(jmax)-corran(1))*365.25+0.5)
!    print*,' BP filt called',jmax,t_small,flow,fhigh,npole,corran(jmax),sig5,mmax,md
!      time0=secnds(0.0)
 
    covar=0.0

    call band_pass_filt(t_small,flow,fhigh,npole,mmax,cr,sig5)
                        
    print*,' BP filt 3',jmax,t_small,flow,fhigh,npole,corran(jmax),sig5,mmax,md
    do  i=1,jmax
      imax=int((corran(i)-corran(1))/t_small+0.5)+1
       sumlast=0.
       do  k=1,i
         dsum=0.0
         kmax=int((corran(k)-corran(1))/t_small+0.5)+1
         do   kk=1,kmax

           dsum=dsum+cr(imax+1-kk)*cr(kmax+1-kk)
         end do
         sumlast=dsum
         covar(i,k)=dsum
       end do
    end do
    print*,' Covariance for BP-filtered noise constructed'


  end subroutine band_pass_cov
      
  subroutine makeMatrix(irow,max,f,covar,nobs,iflip,mmax,nmax)
!  make covariance matrix, covar, from the filter function,f
!  Initially, builds covariance matrix on the assumption unformally sampled observations
!   with no missing data
!  Then, selects rows and columns specified by irow to build actual covariance matrix
!
!     irow  vector contain row numbers (or column numbers) to select out of complete covariance matrxi
!     max  total size of covariance matrix (and f) assuming no missing data
!     f  filter function (length max)
!     covar  covariance matrix
!     nobs  number of data
!     nmax  dimension of covar
!     iflip  if 0; for calculate covariance; if 1--calculating inverse covariance
!     iflip 2  is the CORRECT way of computing the so-called covariance matrix given finv
!
!      Normally, if the H matrix made up of filter coefficents of f (See Langbein 2004, eq 8)
!       Then the cov=H*H'    <------------iflip=0
!     But, if we're given the inverse filter coefficents, finv
!       Then the cov^-1 = H'*H  <---------- iflip=2
    integer, intent(in) :: nmax,mmax,max,nobs,irow(nmax),iflip
    integer :: i,j,k,imax,kmax, nobs2,idif,kcol,kk,kmin,krow
    real(kind=real64), intent(in) :: f(max)
    real(kind=real64), intent(out) :: covar(nmax,nmax)
    real(kind=real64) :: sum,tmp,sum1,sum2
    real(kind=real64), allocatable :: scr(:,:)
    
    allocate(scr(max,max))
!    print*,' max nobs, mmax nmax', max, nobs, mmax, nmax
    if (iflip .lt. 2 ) then
        if (nobs .ne. max ) then
!  do this when there are missing data
          do i=1,max
            imax=i
            do  k=1,i
               sum=0.0
               kmax=k
!         if (( i .ge. 1) .or. (k .ge. 1 )) then
               if (( i .eq. 1) .or. (k .eq. 1 )) then
                 sum1=0.0
                 sum2=0.0
                 kmin=kmax
               else
                 sum=scr(i-1,k-1)
                 kmin=1
               end if
              do  kk=1,kmin

                 sum=sum+f(imax+1-kk)*f(kmax+1-kk)
              end do
              scr(i,k)=sum
              scr(k,i)=sum

            end do
          end do
    
!
!  select row/columns
!
          
          do  i=1,nobs
             krow=irow(i)
             do  j=1,i
               kcol=irow(j)
               covar(i,j)=scr(krow,kcol)
                covar(j,i)=scr(kcol,krow)
              end do
          end do


        else
!  when there are no missing data, do this....  

           do  i=1,max
             imax=i
             do  k=1,i
               sum=0.0
               kmax=k
!         if (( i .ge. 1) .or. (k .ge. 1 )) then
               if (( i .eq. 1) .or. (k .eq. 1 )) then
                 sum1=0.0
                 sum2=0.0
                 kmin=kmax
               else
                  sum=covar(i-1,k-1)
                 kmin=1
              end if
              do  kk=1,kmin

                sum=sum+f(imax+1-kk)*f(kmax+1-kk)
              end do
              covar(i,k)=sum
              covar(k,i)=sum

           end do

         end do   ! index on k
      end if   ! index on i

!
!   flip matrix for inverse covariance
!
!  NOTE, I'm not sure that this is required; tests suggest with iflip=0,
!   that I get an inverse
      if (iflip .eq. 1 ) then
        nobs2=int(nobs/2)
        idif=nobs-2*nobs2
!   idif=0 then, nobs is even; idif=1 then nobs is odd
        do  i=1,nobs2
          do  j=1,nobs
             print*,i,j,nobs+1-i,nobs+1-j
             tmp=covar(nobs+1-i,nobs+1-j)
             covar(nobs+1-i,nobs+1-j)=covar(j,i)
             covar(j,i)=tmp
          end do  ! index on j
        end do ! index on i
        if (idif .ne. 0) then
!          i=nobs2+1
          do  i=1,nobs2+1
          do  j=1,nobs
            print*,i,j,nobs+1-i,nobs+1-j
            tmp=covar(nobs+1-i,nobs+1-j)
            covar(nobs+1-i,nobs+1-j)=covar(j,i)
            covar(j,i)=tmp
          end do    ! index on j
          end do     ! index on i          
        end if

      end if
!  end iflip .ne. 2
      end if
!
!   Do this for inverse covariance matrix given the inverse filter
!
!    NOTE that nobs must equal max
      if (iflip .eq. 2 ) then
        do  i=1,max
!      do 41 i=max,1,-1
          imax=i
          do  k=1,i
!        do 42 k=i,1,-1
            sum=0.0
            kmax=k
!         if (( i .ge. 1) .or. (k .ge. 1 )) then
            if (( i .eq. 1) .or. (k .eq. 1 )) then
               kmin=kmax
            else
               sum=covar(i-1,k-1)
               kmin=1
            end if
            do  kk=1,kmin
!        do 43 kk=kmin,1,-1
!         print*,i,k,kk,imax+1-kk,kmax+1-kk
              sum=sum+f(imax+1-kk)*f(kmax+1-kk)
            end do
            covar(i,k)=sum
            covar(k,i)=sum

          end do

        end do
!
!  flip around the axis from lower left to upper right....
        do  i=1,max
          do j=1,max
            scr(max+1-j,max+1-i)=covar(i,j)
          end do
        end do
        do  i=1,max
          do  j=1,i
            covar(i,j)=scr(i,j)
            covar(j,i)=scr(i,j)
          end do
        end do

   end if
    
 end subroutine makeMatrix


  
end module filterfunc_mod
