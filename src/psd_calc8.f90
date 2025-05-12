program psd_calc8
!  Program computes PSD curve for given input of PSD functions
!    Functions include
!      white noise
!      power law noise
!      power law noise with Gauss-Markov roll-off
!      band pass filtered
!      spiked noise
!  MODIFICATION of psd_cal6 to include computation of PSD
!   using quadrature addition (which is psd_calc) and addition
!   of filter function (as est_noise7.2x using modtype a, f, or c)
!    aka simple addition
!
!  psd_calc8 is a rewrite of psd_cal7 to "modern fortran"


  use iso_fortran_env
  use filterfunc_mod 
  implicit none
!  real(kind=real64), allocatable :: ran,ran
  character(len=1) :: ModType
  real(kind=real32) :: pi
  real(kind=real32) :: tsam,fny,fs,tlen,fl
  real(kind=real32) :: fexp1,amp_pl,amp_pl1,fgm, psd_db, pow_db,P1o
  real(kind=real32) :: fexp2,amp_pl2,P2o,sigma,Pw,alpha
  real(kind=real32) :: flow,fhigh,tmp,amp_bp,h,powBP
  real(kind=real32) :: freq, powPL1, powPL2,power,powSP,pow
  real(kind=real64), allocatable :: ran(:,:),filt(:),wsave(:),dreal(:),dimag(:),tmpx(:)
  real(kind=real64) :: azero
  real(kind=real32) :: sumh
  integer :: npts,i,j,npole,nx,ishift
  pi=2.0*atan(1e+16)
  open(10,file='psdcalc.out')
!  input the stuff needed to create the model PSD
  print*,' Input the sampling interval in days'
  read(5,*)tsam
  fny=0.5*(365.25/tsam)    ! nyquist frequency
  fs=365.25/tsam            ! sampling frequency
  print*,' Input the length of the data set in years'
  read(5,*)tlen
  fl=1/tlen
  npts=int(fny/fl)+1
  
  print*,' Input parameters for first power law noise'
  print*,' Input the index (n=2 for random walk noise)'
  read(5,*)fexp1
  print*,' Input the amplitude'
  read(5,*)amp_pl
  amp_pl1=amp_pl
  if (amp_pl .eq. 0) amp_pl=0.00001
  print*,' Input frequency (radians/yr) for G-M part'
  read(5,*)fgm
  psd_db=10.0*alog10(2.0)-10.*fexp1*alog10(2.0*pi)  &
      +5.0*fexp1*alog10(fs)-10.*alog10(fs)      
  if (abs(amp_pl) .ne. 0) then
     pow_db=20*alog10(amp_pl)+psd_db
     print*,' first power law',pow_db,' db', amp_pl1,fexp1 
     P1o=10**(pow_db/10.)
   else
     P1o=0.0
   end if
   
   print*,' Input parameters for SECOND power law noise'
   print*,' Input the index (n=2 for random walk noise)'
   read(5,*)fexp2
   print*,' Input the amplitude'
   read(5,*)amp_pl
   amp_pl2=amp_pl
   psd_db=10.0*alog10(2.0)-10.*fexp2*alog10(2.0*pi) &
      +5.0*fexp2*alog10(fs)-10.*alog10(fs)  
   if (abs(amp_pl) .ne. 0) then
     pow_db=20*alog10(amp_pl)+psd_db
     print*,' 2nd power law',pow_db,' db', amp_pl,fexp2 
     P2o=10**(pow_db/10.)
    else
      P2o=0.0
    end if      

    print*,' '
    print*,' Input the white-noise standard deviation'
    read(5,*)sigma
    Pw=sigma**2/fny
!    print*,Pw,10.*alog10(Pw)
    
    print*,' '
    print*,' Input parameters for band-pass filtered noise'
    print*,' input the pass band parameters flow, fhigh, c/yr'
    read(5,*)flow, fhigh
    if (fhigh .lt. flow) then
      tmp=flow
      flow=fhigh
      fhigh=tmp
    end if
    print*,flow,fhigh
    print*,' Input the number of poles (1,2,3,or 4)'
    read(5,*)npole
    print*,' Input the amplitude'
    read(5,*)amp_bp

    fgm=fgm/(2.0*pi)   !! convert from radians/yr to cycles/yr
    print*," "
    print*,"  Input the noise type using n, q, or a"
    print*,"    n or q -- quadrature addition"
    print*,"    a -- simple addition"
    read(5,fmt="(a1)")modtype
    
!!!!!!!!!!!!!!!!!!!!!!
!!!!  Program divides here
!!  based upon quadrature or additive noise
!!!!!!!!!!!!!!!     
    if ((modtype .eq. "n" ) .or. (modtype .eq. "q")) then
      print*,"  PSD using quadrature"
      if (amp_bp .ne. 0) then
        sumh=0.0   !!  sum to get variance of BP filtered noise
        do i=1,npts
!  compute variance of bp filtered noise
          freq=fl*float(i)
          h=((freq/flow)**2)    &
           /((1.0+((freq/flow)**2))*(1.0+((freq/fhigh)**2)))
          h=h**npole
          powBP=(4.0*1.0**2)*h
          sumh=sumh+powBP*fl 
        end do
        sumh=1.0/(sumh)  !! normalization to get standard deviation to be 1.0
        print*,' Normalize BP to get 1-unit variance ', sumh
      end if             
      do  i=1,npts
        freq=fl*float(i)
!       power law
        powPL1=P1o
        powPL2=P2o

!  modification per simon williams
       if (fexp1 .ne. 0) powPL1=P1o/((fgm**2 + freq**2)**(fexp1/2.0))
!       if (fexp1 .ne. 0) powPL1=P1o/(fgm**fexp1 + freq**fexp1)
       if (fexp2 .ne. 0) powPL2=P2o/( freq**fexp2)
!  white noise 
        Pw=Pw  
! bandpass filtered
        h=((freq/flow)**2)    &
         /((1.0+((freq/flow)**2))*(1.0+((freq/fhigh)**2)))
        h=h**npole
        powBP=(4.0*amp_bp**2)*h*sumh
          
! spike noise  not used
!        df=abs(freq-fspike)
       powSP=0.
!        if (df .lt. 0.5*fl) powSP=Ps
 
! Sum the noise       
        pow=powPL1+Pw+powBP+powSP+powPL2
        power=10.*alog10(pow)
  
        write(10,*)freq,power
      end do
!      print*,'sumh=',sumh,' fl=',fl,' npts=',npts,' (sqrt(sumh/npts)=',sqrt(sumh/float(npts))
      
      
    else
      print*,"  PSD using additive filters"
!  this is way more complex than the quadrature case
!   first, build the filter in the time-domain, then computes its FFT to output ij frequency domain
!    will need to use subroutines contained in filterfunc_mod

!  force npts be a power of 2
      npts=npts*2
      npts=int(alog(float(npts))/alog(2.))+1  
  
      npts=1*(2**npts)
  
      fl=(365.25/tsam)/float(npts)
      tsam=tsam/365.25
      allocate(ran(4,npts))
      allocate(filt(2*npts))
      allocate(tmpx(npts))
      do j=1,4
      do i=1,npts
        ran(j,i)=0.0d+0
      end do
      end do
!  white noise

      if (sigma .ne. 0.) then
         ran(1,1)=dble(sigma)
      end if
 
!!  power law noise   
    if (amp_pl1 .ne. 0.0) then

       alpha=fgm*(2.0*pi)  !! convert from c/yr to radians/yu
       call frac_diff(tmpx,fexp1,alpha,tsam,npts)
       ran(2,:)=tmpx
       print*,fexp1,alpha,tsam,npts
    
        amp_pl1=amp_pl1*(tsam**(fexp1/4))

        do  i=1,npts
          ran(2,i)=ran(2,i)*amp_pl1
        end do
      end if
      if (amp_pl2 .ne. 0.0) then
        alpha=0.0
        call frac_diff(tmpx,fexp2,alpha,tsam,npts)
        ran(3,:)=tmpx
        amp_pl2=amp_pl2*(tsam**(fexp2/4))
        do  i=1,npts
          ran(3,i)=ran(3,i)*amp_pl2
        end do
      end if
      
!!  BP filtered
      if (amp_bp .ne. 0.) then

        call band_pass_filt(dble(tsam),flow,fhigh,npole,npts, tmpx,amp_bp)
        ran(4,:)=tmpx
      end if
      do i=1,2*npts      !!!!  note that filt dimension is twice the required amount
        filt(i)=0.0d+0
      end do
!!  sum all of the impulse responses and pad first and last half of filt with zeros

      do i=1,npts
        filt(npts/2+i)=0.0d+0
        do j=1,4
           filt(npts/2+i)=filt(npts/2+i)+ran(j,i)
        end do

        if (i .ge. 9*npts/10) then     !!!  cosine taper last 10% of filter
           filt(i+npts/2)=filt(i+npts/2)*cos((pi/2.0)*(i-9*npts/10)/(npts/10))
        end if
!        write(7,*)i,i+npts/2,filt(i+npts/2)
      end do
      
!  start the fft process

      

!  compute fft
      allocate(wsave(3*2*npts+15))
      allocate(dreal(npts))
      allocate(dimag(npts))
      call dzffti(npts*2,wsave)    !! initialize

      call dzfftf(npts*2,filt,azero,dreal,dimag,wsave)   !!!  fft

      do  i=2,(npts)-1
!   do some simple smoothing  (3-points)
        Pow=0.0
        do j=1,3
            Pow=4.0*(dreal(i-2+j)**2 + dimag(i-2+j)**2) + Pow
        end do
        Pow=Pow/3.0
        Pow=10*alog10(Pow/fl) + 10.0*alog10(float(npts)/2.0)
!        Pow=10.0*dlog10(4.0*(dreal(i)**2 + dimag(i)**2)/fl) &
!             + 10.0*alog10(float(npts)/2.)
          write(10,*)fl*float(i)/2.0,Pow
      end do
!      i=(npts)

!      pow=10.0*dlog10((dreal(i-1)**2 + dimag(i)**2)/fl) &
!            + 10.0*alog10(float(npts)/2.)
!      write(10,*)fl*float(i),Pow
    
      deallocate(ran)
      deallocate(filt)
      deallocate(wsave)
      deallocate(dreal)
      deallocate(dimag)
      deallocate(tmpx)
    end if    !!!!  end switch between additive and quadrature models
  close(10)
  print*,' Output PSD model in psdcalc.out'
  print*,'  Freq (c/y), Power, db unit^2/(c/yr)'
end program psd_calc8

