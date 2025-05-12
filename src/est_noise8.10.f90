program est_noise8
!  significant rewrite of est_noise7.30
!   Fortran now is compatible with f90 and more recent version (old versions are fortran 77)
!  Main motivation for the re-write is memory allocations; older versions with 'largish'
!    arrays did not compile and run.
! New version will be using allocatable arrays.
!  Version 8.1  May 2025
!  1)  found small bug in the frac_dif code .... a few variables needed to be 64-bit reals and not 32-bit
!  2)  Introduced openmp paralleliation on many do-loops. In particular, this increased
!       the cpu speed in the convolution subroutine which is heavily used in the so-called 'f'/fast mode.

!  Notes from f77 version
!C  Revise est_noise6.50 and est_noise7.10 to combine both Max Likelihood estimator (MLE)
!c   to determine noise models to time-series data and Simultaneously, fit a simple temporal
!c  function to the time-series.
!c
!c  Noise model for est_noise6.50 adds several independent noise sources together as
!c
!c    N(t)^2 = (h_1 * e_1)^2 + (h_2 * e_2)^2 + (h_3 * e_3)^2 ...
!c
!c   where * is convolution, h_i are impluse responce of a i_th noise model and e_i are
!c     independent white noise sources
!c
!c  An alternate noise model consists of simply adding the impluse responses and then convolving with
!c   a single white noise source;
!c
!c   N(t)=  [ h_1 + h_2 + h_3 +...] * e(t)
!c
!c   which is used in est_noise7.10 and its variants
!c
!c   The additive form as the advantage of rapid compution of the inverse of the data
!c   covariance, and for missing data, using the method of Bos et al (2012) to speed up
!c   the inversion.
!c
!c   Currently, noise models are
!c   Power-law noise as generated using Hosking (1981)
!c   Gauss Markove noise discussed by Langbein (2004)
!c   white noise
!c   Band Pass filtered noise, Langbein (2004)
!c   Program could be expanded to ARIMA models using Hosking (1981)
!c
!c  Temporal dependence includes
!c   Nominal value including one or more time-series
!c   rate
!c   changes in rate over a specified time
!c   offsets at prescribed times
!c   sinusoids (in phase and quadrature terms)
!c   exponential trends
!c   Omori law (log(1+t/tau)
!c   User prescibed external functions (Pressure data for strainmeter adjustment)
!!c  Program is configured such that standard input used historically for est_noise6.50 will
!c   work as default for this new version.
!c
!c  Program includes option to create noisy data sampled at time of real data for simulations
!c
!c  Program includes additional time-stamps other than those used befoer;
!!c   year day-of-year
!c   year, month, day (numbers)
!c   GMT  (Generi!c Mappping Tools time stamps)
!c   MJD (Modified Julian Day)
!c
!c  Command to compile on Ma!c OS 10.6 with Intel compiler and MKL libraries
!c  ifort   -Wl,-stack_size,0x40000000 -O3      -L/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal -lmkl_lapack -lmkl_core -lmkl_intel_lp64 -lmkl_intel_ilp64 -lmkl_intel_thread -o est_noise7.22 est_noise7.22.f ~/proglib/sublib/time.f funmin7.22.f invfft.f NedlerMeadSimplex.f  modify.f NedlerMeadSimplexEXP7.22.f
!c
!c  Using gfortran
!c  gfortran -O3 -o est_noise7.21 est_noise7.21.f time.f  funmin7.21.f invfft.f  -framework Accelerate NedlerMeadSimplex.f  modify.f NedlerMeadSimplexEXP7.21.f
!c
!c  Modify Apr 2021 to correctly compute covariance matriculates when encountering GGM noise model
!c
!  revised stuff for version 8 (August/November 2024)
! 1)  Program can accomodate dates between 1800 and 2200 rather than 1960 to 2040
! 2) fixed error bar calculation for the noise parameters
! 3)  New normalization for BP filter -- sums of squares of filter coefficients now equal 1.0
! 4) Program suite now updated to more modern Fortran which provides active memory management
!     through allocate function. Version 7 stopped compiling and running using more modern compilers.
!      This caused a significant re-write of the programs.

! DECLARATUION
  use iso_fortran_env
  use alloc_mod
  use est_noise_mod
  use GetStuff_mod
  use genNoise_mod
  use time_mod
  use filterfunc_mod
  use mle_mod
  use NedlerMeadSimplex_mod
  use model_fit_A_mod
  use omp_lib


  implicit none

  character(len=3) :: netx(20),netp                    !!  date format
  character(len=1) :: ans, ModType                      !! y/n
  character(len=132) :: filename, string
  character(len=6) :: nameNoise(7)
!  integer :: nbt, nrate, n_rat_chg, nper, noff, n_exp, 
  integer:: nobs(10)
!  integer :: nexp_fix,n_exp_fix, n_file_press,nc   !!  number of model parameters
  integer :: i,j,k,kz, nopt,nc                                 !!  indices for do-loops
!  real(kind=real32) :: rate_norm                !!  normalization of rate
!  real(kind=real64) ::  t_start,t_stop          !!  start and stop times in num of days since 1/1/2800
!  real(kind=real64) :: rat_chg1(:),rat_chg2(:), rat_chng_norm(:)
!  real(kind=real64) ::  per(:),off(:),texp(:),aux_norm(:)      !! model descriptions
!  real(kind=real32)  ::  bexp(:)
  character(len=7) :: choice(10)
  integer ::  max_parm=82,  ios, iln,loop, iFlagReLoop,loopOrig,iloop,iswt,nx,icx

  logical :: filexist  
  real(kind=real32) :: dist,err,tsmall,t_lag,press_min,press_max
  real(kind=real32) :: ambp(8,7),amby(8),dith(7)
  real(kind=real32) ::fmiss,fpress,dt_ave,dt_max,dt_min, dt_sd,Tstart
  real(kind=real64) :: timex,pi,txxx,chi2,sum_press
  real(kind=real64) :: tpress,ZBQLNOR
  real(kind=real64), allocatable ::  tnoise(:)      ! used if random number generator is used     
  real(kind=real32), allocatable :: dnoise(:)        ! used if random number generator is used  
  real(kind=real64) :: dec_time,fjul,fmle,fmax,fv1,fv2,fv3,fv4
  integer :: itlen,jul,itime,nyr,ierr,nchoice,ii,jj,kkk
  real(kind=real32) :: est_white_noise,sum_diff,perNmissEst,det,f_mle,rms,sig1,sig_color,siginstr,dif,fmin,fxx
  real(kind=real32) :: alpha,alpha2,amp1,amp2,amp_bp,eexp,exp1,exp2,f1a,f2a,flow,fhigh,PPo,xamp,wh_add,xfl,xfny,timex0,timex1
  integer :: n_diff,idec,nmissEst,nptsThr,idiff,ix,ic_orig,IOstatus,npole, &
       maxk,ifloat(8),nfloat,ITER,irow,isave,count0, count_rate, count_max,count1
  real(kind=real64), allocatable :: A_orig(:,:), time_orig(:), t_year_orig(:)
  real(kind=real32), allocatable :: d_orig(:)
  real(kind=real32) :: amp_mod,damp1,dsig,freq,fs,h,powBP,psd_db,psd_hf,psd_hf_pl1, &
     psd_hf_pl2,psd_pl,psd_pl1,psd_pl2,psd_wh,sig_min,alpha_2pi,dalpha,damp_bp,dexp1,dexp2,flong,damp2,FTOL,small, &
     xn(7),xx(7),enull
  real(kind=real32) vsig(8),vamp1(8),vexp1(8),vamp_bp(8),valpha(8),vamp2(8), vexp2(8),acov(7,7),cov(7,7),wz(23)
  integer :: numThread,maxThread
  pi=2.0d+0*datan(1d+32)
  data (vsig(i),i=1,8)/0.7,1.2,1.2,1.2,0.7,1.2,1.2,0.7/
  data (vamp1(i),i=1,8)/1.95,1.95,0.3,0.3,0.3,0.3,1.95,1.95/
  data (vexp1(i),i=1,8)/0.4,-0.4,-0.4,-0.4,-0.4,-0.4,0.4,-0.4/
  data (valpha(i),i=1,8)/4.5,0.35,4.5,4.5,4.5,4.5,0.35,4.5/
  data (vamp_bp(i),i=1,8)/0.55,1.79,0.55,1.79,0.55,1.79,0.55,0.55/
  data (vamp2(i),i=1,8)/2.31,0.26,2.31,0.26,0.26,0.26,2.25,0.26/
  data (vexp2(i),i=1,8)/-0.6,-0.6,0.6,0.6,-0.6,-0.6,0.6,0.6/
  data (nameNoise(i),i=1,7)/'white ','plamp1','plexp1','GM    ','BPamp ','plamp2','plexp2'/
  
   Anor=1.0
!  start wall clock and cpu clock
  call system_clock(count0, count_rate, count_max)
  call cpu_time(timex0)

!  get number of threads (openmp) available and limit
!  that number to 8
   maxThread=OMP_get_max_threads()
   print*,' Max number of threads available',maxThread
   if (maxThread .gt. 8 ) then
     call OMP_set_num_threads(8)
     maxThread=OMP_get_max_threads()
     print*,' change thread count to', maxThread
   end if
!!  check to see if there is a seed file, seed.dat
 

  open(23,file="seed.dat")
  read(23,*,iostat=ios)iseed
  if ( ios .ne. 0) then
    print*, 'program will use its default value'
    print*, ' best to create a seed file before execution'
    print*, 'echo positive_integer > seed.dat'
    iseed=728745902
    rewind(23)
    write(23,*)iseed

  end if

 
!!  Open statements
  open(25,file='estin.jrn',action='write')   !!  Journal file of input to est_noise
  open(89,file='tauexp.out',action='write') !! file showing various trial of estimating the 
!!                                           time-constants of the data if time constants are "float"
  open(17,file='adj.out',action='write')      !!  output can be used to drive the adjust program
  open(12,file='max.out',action='readwrite')      !!  output of iterating to maximize MLE
  open(13,file='model.out',action='readwrite')      !!  LS solutions of trajectory model for iterating to maximize MLE
!  open(12,file='max.out')      !!  output of iterating to maximize MLE
!  open(13,file='model.out')      !!  LS solutions of trajectory model for iterating to maximize MLE
  open(14,file='model.dat',action='write') 
  open(15,file='max.dat',action='write')  
  print*,' A journal file consisting of typed input'
  print*,'  is created in estin.jrn'
  
!!  Ask lots of questions about data and the analysis that follows
  print*,' Program estimates the power spectral density (PSD)'
  print*,'  of data.  PSD functions can be a combination of:'
  print*,'  1) white noise'
  print*,'   2) Simple power law noise P/f^n'
  print*,'  3) a Gauss Markov version P/(fa^n + f^2)'
  print*,'      where fa=alpha/2*pi'
  print*,'  4) A second power law'
  print*,'  5) Band-pass filtered noise'
  print*,' '
  print*,' While estimating the PSD function, program also estimates'
  print*,'  various parameters that describe the time series including'
  print*,'  1) DC term'
  print*,'  2) rate'
  print*,'  3) sinusoidal amplitudes of specified frequencies'
  print*,'  4) rate changes'
  print*,'  5) offset'
  print*,'  6) simple exponential and log(t) function'
  print*,'  7) user supplied function or data'
  print*,' '
  print*,' Program can handle data in various formats'
  print*,' Version 8.10  --- May, 2025'
  print*,' '

  print*,' A journal file consisting of typed input'
  print*,'  is created in estin.jrn'

  print*,'Input the data type for  processing'

  print*,' doy=data  or otr=data'
  print*,'       format of year, day of year, data, error bar'
  print*,' otd=data'
  print*,'       format of YearMnDa,  data, error bar'
  print*,' otx=data'
  print*,'       format of year mo da,  data, error bar'
  print*,' mjd=data'
  print*,'       Modified Julian day, data, error bar'
  print*,'For all of the above, the da is double precision'
  print*,'    which allows decimal days to accomodate seconds'
  print*,' '
  print*,' gmt=data with GMT format for time'
  print*,'        year-mo-daThr:mn:secs.x data, error bar'
  read(5,*)net
  if ((net .ne. 'mjd')  .and.  (net .ne. 'otr') .and. (net .ne. 'otd') .and. (net .ne. 'gmt') &
      .and. (net .ne. 'otx') .and. (net .ne. 'doy')) then
        print*,'Must be either doy otr otx otd mjd or gmt'
        stop
  end if
  write(25,*)net
!
!   Read Inputs that specifies the time-dependent model
!
  nmod=0
  print*,' input number of time series'
  print*,'   program will estimate only one set of PSD functions'
  read(5,*)nbt
  write(25,fmt="(i4, '    # number of time series')")nbt
  nmod=nmod+nbt
  print*,' input the period of interest'

  call GetTime1(net,t_start,t_stop,5,1)
  print*,' Day numbers from 1800 ',t_start,t_stop
  print*,' '
  print*,' Input the parameters of time series to be estimated'
  print*,' '
  print*,' Will rate be estimated? y/n' 
  read(5,*)ans
  write(25,*)ans
  nrate=0
  if (ans .eq. 'y') then
    nrate=1
    rate_norm=((t_stop-t_start)/2.)/365.25
    print*,' rate renomalization is ',rate_norm
  end if
  nmod=nmod+nrate
  print*," "
  print*,' Rate changes:'
  print*,'  time is the "hinge point"'
  print*,' Input the number of rate changes in data'
  read(5,*)n_rat_chg
  write(25,fmt="(i4,'    #  Number of rate changes')")n_rat_chg
  nmod=nmod+n_rat_chg
  if (nmod .gt. max_parm) then
     print*,' Number of parmeters exceed limit of', max_parm
     stop
  end if
  if (n_rat_chg .ne. 0) then
    do i=1,n_rat_chg
      call GetTime2(net,i,rat_chg1(i), rat_chg2(i),5,1)
      if (rat_chg1(i) .ge. t_stop) then
         print*,' rate change beyond end of time series'
         stop
       end if
       rat_chng_norm(i)=(rat_chg2(i)-rat_chg1(i))/365.25
       print*,' Interval of rate change', i
       print*,rat_chg1(i), rat_chg2(i)
       print*,' Normalization of rate change is:', rat_chng_norm(i)
    end do
 end if
 print*,' '
 print*,' Input the number of periodicities in the data'
 read(5,*)nper
 write(25,fmt="(i4,'    # Number of periodicities')")nper
 nmod=nmod+nper*2
 if (nmod .gt. max_parm) then
    print*,' Number of parmeters exceed limit of', max_parm
    stop
 end if
 if (nper .ne. 0) then
    do  i=1,nper
      print*,' Input the period (in days) of ', i,' period'
      read(5,*)per(i)
      write(25,*)per(i)
    end do
 end if
 print*,' '
 print*,' Offset:  specify first day after offset'
 print*,' Input the number of offsets'
 read(5,*)noff
 write(25,fmt="(i4,'    # Number of offsets')")noff
 nmod=nmod+noff
 if (noff .gt. max_parm) then
    print*,' Number of parmeters exceed limit of', max_parm
    stop
 end if
 if (noff .ne. 0) then
   do i=1,noff
     call GetTime3(net,i,off(i),5,1)
     print*,' Offset at ',off(i)
     if (off(i) .gt. t_stop) then
        print*,' Offset beyond end of time series'
        stop
     end if
   end do
 end if
 print*,' '
 print*,' Input the number of exponentials in time series'
 read(5,*) n_exp
 write(25,fmt="(i4,'    # Number of exponentials/log func.')")n_exp

 nmod=nmod+n_exp
 nexp_fix=0
if (n_exp .ne. 0) then
      do  i=1,n_exp
        call GetTime4(net,i,texp(i),5,1)
!        texp(i)=texp(i)-t_start
        print*,' Exponential starting at ',texp(i)
        
        if (texp(i) .ge. t_stop) then
           print*,' exponential beyond end of time series'
           stop
        end if     
        print*,' Input the time constant in years and fix/float'
        read(5,*)bexp(i),exp_choice(i)
        write(25,*) bexp(i),exp_choice(i)
        print*,' Input the type of function; e for 1-exp(-t/tau)'
        print*,'      or                     m for log10(1.0 + t/tau)'
        read(5,*)exp_type(i)
        write(25,*)exp_type(i)
      end do
      nexp_fix=0
      do  k=1,n_exp
        if (exp_choice(k) .eq. 'fix    ') nexp_fix=nexp_fix+1
      end do
    print*,' The number of exponential time constants that are fixed:', nexp_fix
 end if
   
!
!   Input the time series.  This happens in two stages; first-- query each data file for
!     the number of observations, then second, read the data and construct A matix
!   To keep the inputs between older versions of est_noise, the two steps are
!     separated with addition input concerning auxilary data (ie, atmospheric pressure)
!
!   Part one of data read
 ic=0
 do i=1,nbt    !!! loop over number of time-series (ie, baselines)
    print*,' Input the format style of time-series data (otr, doy, otx, otd, mjd, or gmt)'
    read(5,*)netx(i)
    write(25,*)netx(i)
    print*,' Input name of file for time-series number ',i
    read(5,fmt="(a132)")filename
    write(25,fmt="(a132)")filename
    inquire(file=filename,exist=filexist)
    if ( filexist  ) then
       print*,' data file is present'
    else
       print*,' data file ', filename, ' does not exist'
       print*,'  Balling'
       stop
    end if   !! fileexist?
    open (20+i,file=filename,action='read')
 !!  get number of data
    do 

      read(20+i,fmt="(a132)",end=10) string
      call GetData(netx(i),string,timex,dist,err)
      if ((timex .ge. t_start) .and. (timex .le. t_stop)) ic=ic+1
    end do  !! end short loop to read and count data
10  continue  
    print*,' number of observations',ic 
!    close(20+i)
 end do     !!  finish looping over nbt num of time series

!  Now, input the auxiliary (pressure data); will return to data read
  print*,' Number of files of Auxillary data (pressure)'
  read(5,*)n_file_press
  write(25,fmt="(i4,'    # Number of auxillary time series')")n_file_press
!
!  Allocate sizes of various matrices
!

!  allocate(A(ic,nmod+n_file_press))
!  allocate(t(ic))
!  allocate(d(ic))
!  allocate(res(ic))
  nmod=nmod+n_file_press
  call CreateGlobal
  nmod=nmod-n_file_press
!!
!  Rewind the data file; read the data and populate the design/trajectory matrix A
!
  ic=0
  open(26,file='prob1.out')
  do i=1,nbt   !!! loop on each time-series file
     rewind(20+i)
     do 
       nobs(i)=0
       read(20+i,fmt="(a132)",end=11) string
       call GetData(netx(i),string,timex,dist,err)     !!!  note -- err is ignored
       if ((timex .ge. t_start) .and. (timex .le. t_stop)) then
         ic=ic+1
         nobs(i)=nobs(i)+1
         t(ic)=timex-t_start    !!!   reference time is now t_start
         d(ic)=dist
!  nominal value
         nc=0
         A(ic,i)=1.0d+0
         nc=nc+nbt
!  rate
         if (nrate .eq. 1) A(ic,nc+1)=(t(ic)/365.25)/rate_norm
         nc=nc+nrate  
!  rate changes
         if (n_rat_chg .gt. 0) then
           do  k=1,n_rat_chg
             A(ic,nc+k)=0.0
             if ((timex .ge. rat_chg1(k)) .and.  (timex .lt. rat_chg2(k))) &
               A(ic,nc+k)=((timex-rat_chg1(k))/365.25)/rat_chng_norm(k)
    
             if (timex .ge. rat_chg2(k))  &
                 A(ic,nc+k)= ((rat_chg2(k)-rat_chg1(k))/365.25)/rat_chng_norm(k)
           end do  
         end if    !! end rate change
        nc=nc+n_rat_chg  
!! periodicities
        if (nper .gt. 0) then
          do  k=1,nper

            txxx=2.0d+0*pi*t(ic)/per(k)
            A(ic,nc+2*k-1)=cos(txxx)
            A(ic,nc+2*k-0)=sin(txxx)
          end do
        end if
        nc=nc+2*nper
!  offsets in data
        if (noff .gt. 0) then
          do  k=1,noff
            A(ic,nc+k)=0.0
            if (timex .ge. off(k)) A(ic,nc+k)=1.d+0
          end do
        end if
        nc=nc+noff
!
!  Exponentials
!
 
       if (n_exp .ne. 0 ) then
       kz=0
         do  k=1,n_exp
           if ( exp_choice(k) .eq. 'fix') then
             kz=kz+1
             A(ic,nc+kz)=0.0
             if (timex .gt. texp(k)) then
               if (exp_type(k) .eq. "e") A(ic,nc+kz)=1.0-exp(-(timex-texp(k))/(bexp(k)*365.25))

               if (exp_type(k) .eq. "o") A(ic,nc+kz)=alog10(abs(bexp(k)) +  &
                                           (sngl(timex-texp(k))/365.25))
               if (exp_type(k) .eq. "m") A(ic,nc+kz)=alog10(1.0 +           &            
                                   (sngl(timex-texp(k))/365.25)/abs(bexp(k)) )
             end if
           end if
         end do
         
         nexp_fix=kz
       end if  ! end exp      
       nc=nc+nexp_fix  
       write(26,fmt="(i5,f15.5,f10.2,40(3x,e13.4))")ic,t(ic),d(ic),(A(ic,j),j=1,nc)
       end if    !!!  finish populating the data and A matrix
       
     end do 
11 continue    !!  end of data file
   close(20+i)
   end do    !!! end loop on each time series file
   close(26)
   
!!!  input the auxiliary data if requested
  print*,' column of A matrix prior to input pressure data', nc
  print*,' Number of files of Auxillary data (pressure)'
!  compute the tolerance in time to match the pressure and distance data

  tsmall=(t(ic)-t(1))/(float(ic-1)*10.)
!  tsmall=((A(ic,nbt+1) - A(1,nbt+1)))*365.25*rate_norm/(float(ic-1)*10.)
  if (n_file_press .gt. 0) then

     open(88,file='dist_press.dat',action='write')    !!!  file that contains matching pressure and "distance" data
     print*,'size A',size(A,dim=1),size(A,dim=2),size(d),size(t)
     do i=1,n_file_press

       print*,' Input the format style of baseline data',' (otr,doy, gmt, mjd, otx, otd)'
       read(5,*)netp
       write(25,*)netp
       print*,' Input name of file for aux data number ',i
       read(5,fmt="(a132)")filename
       write(25,*)filename
       inquire(file=filename,exist=filexist)
       if ( filexist  ) then
         print*,' pressure data file is present'
       else
          print*,' pressure file ', filename, ' does not exist'
          print*,'  Balling'
       stop
       end if   !! fileexist?
       open (22,file=filename)
       print*,' Input the lag in days'
       read(5,*)t_lag
       write(25,fmt="(e12.5, '    #  lag time in days')")t_lag  
       press_min=999999999.
       press_max=-999999999.
       sum_press=0.0d+0
       fmiss=-9999.0e+20
!  put missing data symbol in pressure column of A matrix
       do  k=1,ic
         A(k,nc+i)=fmiss
       end do
       k=0
       do
         read(22,fmt="(a132)",end=48)string
         
         call GetData(netp,string,timex,fpress,err) 
         tpress=timex-dble(t_lag)-t_start
         
         if ((timex-dble(t_lag) .ge. t_start) .and. (timex-dble(t_lag) .le. t_stop)) then        
        
           do j=1,ic
!   scan through the distance data to match times of measurement
             if (abs(tpress - t(j)) .lt. dble(tsmall)) then
                
                A(j,nc+i)=fpress
 !               write(88,fmt="(i5,2f15.10,2x,2f10.3)")j,tpress,t(j),fpress,d(j)
                write(88,fmt="(i5,2f15.10,2x,2f10.3)")j,timex-dble(t_lag)-t_start,t(j),fpress,d(j)

                if (fpress .gt. press_max) press_max=fpress
                if (fpress .lt. press_min) press_min=fpress
                sum_press=sum_press+fpress
                k=k+1
             end if
          end do

          
         end if   !!! is time of press data fall between t_start and t_stop
       end do
48     continue
       
       print*,' number of press readings',k
       print*,' min. max, and average pressure',press_min,press_max,sum_press/float(k)
!!  normalize pressure data and replace any missing (fmiss) elements with 0.0
       aux_norm(i)=dble((press_max-press_min))
       print*,' Renormalize auxilary data set ',i,' by',aux_norm(i)
       do j=1,ic
         if (A(j,nc+i) .eq. fmiss) then
            A(j,nc+i)=0.0d+0
         else

            A(j,nc+i)=(A(j,nc+i) - sum_press/float(k))/aux_norm(i)
         end if
       end do
         
  
     end do   !!  end loop for pressure files  
     close(88)
  end if    !! end inputting pressure data
   
   
!  end input of auxillary data

      nmod=nbt+nrate+n_rat_chg+2*nper+noff+n_file_press+nexp_fix
      nmod_orig=nmod
     print*,' Number of data read is ',ic
 !     print*,' Number of model parameters is ',nmod
 !     write(6,*)' Number of model parameters is ',nmod
 !   print*,nmod

!  put data into chonological order
  call chron
!
  open(7,file='prob2.out')
  do  i=1,ic
    write(7,fmt="(i5,f15.5,f10.2,40(3x,e13.4))")i,t(i),d(i),(A(i,j),j=1,nmod)
  end do
  close (7)
  print*, ' '
  print*,'  prob1.out; time, data, and A matrix;  aux data excluded'

  print*,'  prob2.out includes aux data and data sorted into chronological order'
      

!  get statistics on sampling interval

  dt_ave=(t(ic)-t(1))/float(ic-1)
  dt_min=999999.
  dt_max=-99999.
  dt_sd=0
  do  i=2,ic
   if (t(i)-t(i-1) .gt. dt_max) dt_max=t(i)-t(i-1)
   if (t(i)-t(i-1) .lt. dt_min) dt_min=t(i)-t(i-1)
   dt_sd=dt_sd+ (t(i)-t(i-1) -dt_ave)**2
  end do
  dt_sd=sqrt(dt_sd/float(ic-2))
  print*, ' '
  print*,' Average Sampling interval:  ',dt_ave, ' days'
  print*,' Standard deviation of Ave. interval: ',dt_sd,' days'
  print*,' Shortest Sampling interval: ',dt_min, ' days'
  print*,' Longest Sampling Interval:  ',dt_max, ' days'
  print*,' '
  print*,' Input the minimum sampling interval in days to use'
!  read(5,fmt="(d20.15)")dt_sam
  read(5,*)dt_sam
  write(25,fmt="(f20.15,  '    # sampling inteval in days') ")dt_sam
  
  if (dt_min .lt. 0.1*dt_ave) then
      print*," The following may be near simultaneous observations"
      print*,"   and on might need be eliminated"
      do  i=2,ic
        if ( t(i)-t(i-1) .le. 0.1*dt_ave) then
          write(6,5302)i-1,i,t(i-1),d(i-1),t(i),d(i),t(i)-t(i-1) 
        end if
      end do
      stop
  end if
5302  format("i-1 and i",2i5," t(i-1),d(i-1)",f18.9,f9.2," t(i),d(i)",f18.9,f9.2,"  dt=",f12.9)



  print*,' '
  print*,' Input the type Noise model construction'
  print*,' n or q = noise components are summed as squares'
  print*,'   as in est_noise6.x (Langbein, 2004)'
  print*,' a = noise component are simple sums'
  print*,'   automatically switches between cholesky decomposition'
  print*,'    or computing inverse using combo of inverse convoluton and Bos et al (2012) for missing data'
  print*,' c = noise component simple sums'
  print*,'     forces inversion using cholesky decomposition'
  print*,' f = noise component simple sums'
  print*,'    forces combo of inverse convolution and and Bos et al (2012) for missing data'
  print*,' '
  read(5,*)ModType
  write(25,*)ModType
  if ( ModType .eq. 'q') ModType='n'
  print*,' Do you want to substitute real data with random numbers? y/n'
  read(5,*)ans
  write(25,*)ans
  if (ans .eq. 'y') then     
!
!  substitute noise for real data; noisy data has same sampling interval as real data
!   NOTE -- This only works correctly if a single baseline is being examined! (nbt=1)
!
    open(77,file='noise.dat')
!    call genNoise(ModType,t,d,ic,max_data,iseed,dt_sam,max_time)
!  compute length of time series would be if there were no data gaps
    itlen=int((t(ic)- t(1))/dt_sam + 1) + 1
    print*,' number of noisy data without gaps; program will create a set of data at time of actual measurements'
    print*,itlen
    allocate(tnoise(itlen))
    allocate(dnoise(itlen))

    call genNoise(ModType,tnoise,dnoise,iseed,dt_sam)
!  re-sample the continuous data vector -- revised data are those corresponding to time of actual measurements
    do  i=1,ic
      j=int((t(i))/dt_sam + 1)
    
      itime=int(t(i)+t_start)
      dec_time=(t(i)+t_start)-int(t(i)+t_start)
      call num2doy(itime,nyr,jul)
      fjul=float(jul)+dec_time
      d(i)=dnoise(j)
      write(77,*)nyr,fjul,d(i),t(i)/365.25
 
    end do
    deallocate(tnoise)
    deallocate(dnoise)
    print*,' random noise data written to noise.dat'
    close(77)
  end if   !! finish creating random noise for data

!!  estimate the white noise component by taking the 1st differences between data 
!!  and save "original A, t, and d" 
  
  allocate(A_orig(ic,maxnmod))
  allocate(time_orig(ic))   
  allocate(d_orig(ic))
  allocate(t_year_orig(ic))

  ic_orig=ic
  sum_diff=0.
  n_diff=0
  do i=1,ic
     t_year(i)=(t(i)-t(1))/365.25d+00
     t_year_orig(i)=t_year(i)
     d_orig(i)=d(i)
     time_orig(i)=t(i)
     if (i .ne. 1) then
        if (t(i) - t(i-1) .le. dt_sam*1.01) then
           sum_diff=sum_diff+(d(i)-d(i-1))**2
           n_diff=n_diff+1
        end if
     end if
     do j=1,nmod
      A_orig(i,j)=A(i,j)
    end do
  end do


  est_white_noise=sqrt(sum_diff/float(n_diff))/sqrt(2.0)
  print*,' Estimate of white noise component of data is:', est_white_noise

  
!
!  decimate the data!  ---  IMPLIMENTED  only for ModType=n 
!
  print*,' methods to decimate data'
  print*,' 0 = no decimation'
  print*,' 1 = keep 2, skip 1, keep 2, skip 1'
  print*,' 2 = keep 2, skip 1, keep 1, skip 2, keep 2, skip 1...'
  print*,' 3 = keep 2, skip 1, keep 1, skip 2, keep 1, skip 3,',' keep 2 skip 1...'
  print*,' 4 or more...make it option 3'
  print*,'input choice, BUT if the noise model is simple sums,', ' then ignore'
  read(5,*)idec
  if (ModType .ne. 'n') idec=0       !! only implimented for ModType =n
  if (idec .ge. 4) idec=3
  write(25,fmt="(i4,'     # Decimation type')")idec

  if (idec .ne. 0 ) then
     call decimate(idec)
     open (77,file='decimate.dat')
      do i=1,ic
        itime=int(t(i)+t_start)
        dec_time=(t(i)+t_start)-int(t(i)+t_start)
        call num2doy(itime,nyr,jul)
        fjul=float(jul)+dec_time
        write(77,*)nyr,fjul,d(i),t(i)/365.25
      end do
      close(77)
      print*,' Decimated data in decimate.dat'        

   end if
   
!
!  Estimate the data completeness
!
  nptsThr=int((t(ic)-t(1))/dt_sam +0.005)
  nmissEst=nptsThr-ic+1
  perNmissEst=100.*(float(nmissEst)/float(nptsThr))

  print*,' Estimated number of missing data is ', nmissEst
  print*,' Estimated percenta of data gaps is', perNmissEst
  if ( ModType .eq. 'a' ) then
     ModType='f'
!     if (perNmissEst .gt. 25.0) ModType='c'
     if (perNmissEst .gt. 300.*(float(nptsThr))**(-0.25)) ModType='c'
!  empirical formula to estimate whether the fastest method uses 'f' or 'c'
      print*,' ModType changed from a-automatic, to ',ModType
  end if



!  call CreateCovMatrix(ic,ic)
  do i=1,ic
    covar(i,i)=1.0d+0
    covinv(i,i)=1.0d+0
  end do
  nopt=1
  
  call CreateModTypeQ
  call model_fit(nopt,'n',chi2)
  call DestructModtypeQ
  rms=sqrt(chi2/float(ic-nmod))
  print*,' RMS fit using white noise model is',rms
! compute MLE for the white-noise assumption for covariance
  det=float(ic-0)*alog(rms)
  f_mle=-det-0.5*chi2/(rms**2)-0.5*float(ic)*alog(2.0*3.14159265)
  print*,' Log MLE for white noise error model is ',f_mle
  print*, ' '
!  modify white noise version of error bars by scaling by rms
  do i=1,nmod
     ex(i)=ex(i)*rms
  end do

!  Output the results from the  white noise mode
  call output


!  increase dimension of A, t, d,t_year to span all time, even for times 
!   missing data 
   itlen=int((t(ic)- t(1))/dt_sam + 1) + 1
   print*,ic,size(d),size(t)
   print*,'itlen=',itlen,ic,t(ic),t(1),dt_sam
   deallocate(A)
   deallocate(t)
   deallocate(d)
   deallocate(t_year)
   deallocate(covar)
   deallocate(covinv)
   deallocate(res)

 
!   allocate(A(int(itlen*1.00),nmod))
   allocate(A(int(itlen*1.00),maxnmod))
   allocate(t(int(itlen*1.00)))
   allocate(d(int(itlen*1.00)))
   allocate(t_year(int(itlen*1.00)))
   allocate(covar(int(itlen*1.00),int(itlen*1.00)))
   allocate(covinv(int(itlen*1.00),int(itlen*1.00)))
   allocate(res(int(itlen*1.00)))
   do i=1,ic
     do j=1,nmod
       A(i,j)=A_orig(i,j)
     end do
     t(i)=time_orig(i)
     d(i)=d_orig(i)
     t_year(i)=t_year_orig(i)
   end do
  
!
!  For ModType  simple sums (ModType .ne. n)
!  Create a revised design matrix (and data) that fills in data-gaps.
!     force rows of revised design matrix to be filled with zeros
!

 nmiss=0

 if (ModType .eq. 'f') then

 
    do  i=1,ic
      if (i .gt. 1) then
         idiff=int(((time_orig(i)-time_orig(i-1))/dt_sam)+0.005)
          if (idiff .gt. 1 ) then
!            print*,"i,idiff,dtime(i),dtime(i-1)",i,
!     &       idiff,dtime_orig(i),dtime_orig(i-1)

            do k=1,idiff-1
               nmiss=nmiss+1
               ix=int(((time_orig(i)-time_orig(1))/dt_sam)+0.005)+1
!             irowmiss(nmiss)=ix-(idiff-1)+k-1+nmiss0
               irowmiss(nmiss)=ix-(idiff-1)+k-1
!              print*,"i,ix,nmiss,irowmiss",i,ix,nmiss,irowmiss(nmiss)
               d(i+nmiss-1)=0.
               t(i+nmiss-1)=time_orig(i-1)+float(k)*dt_sam
               t_year(i+nmiss-1)=t_year_orig(i-1)+float(k)*dt_sam/365.25d+0
               do j=1,nmod
                 A(i+nmiss-1,j)=0.
               end do
             end do
             d(i+nmiss)=d_orig(i)
             t(i+nmiss)=time_orig(i)
             t_year(i+nmiss)=t_year_orig(i)
             do j=1,nmod
               A(i+nmiss,j)=A_orig(i,j) 
             end do      
           else
        
             d(i+nmiss)=d_orig(i)
             t(i+nmiss)=time_orig(i)
             t_year(i+nmiss)=t_year_orig(i)
             do  j=1,nmod
               A(i+nmiss,j)=A_orig(i,j) 
             end do     
           end if
         end if
    end do
    print*,' Number of missing data are', nmiss
    ic=ic+nmiss
    print*,' percentage of data gaps', 100*(float(nmiss)/float(ic))

        
  end if   !!  branch on 'f'
  print*,' sampling interval in yrs ',dt_sam/365.25
  print*,' '  
  
  
!!  Input the guess for a noise model
  sig1=rms
  print*,'  Input the initial parameters of the PSD and whether the item is "fix" or "float"'
  print*,'   if "fix", then the item is not estimated'
  print*,'   if "float", then the item is estimated'
  print*,' '
  print*,' Input the white noise "instrument precision and fix/float'
  print*,'  Suggested input is ', est_white_noise,' enter -99999.0 to use suggested value'
  read(5,*)siginstr, choice(1)
  write(25,fmt="(f8.3,1x,a7,'    # white noise')") siginstr, choice(1)
  if (siginstr .eq. -99999.0 ) print*,' will try ',est_white_noise,' as white noise component'
  if (siginstr .eq. -99999.0 ) siginstr=est_white_noise
  print*,' Input the amplitude first Power law function and fix/float'
  print*,' To use default value, enter -99999.0'
  read(5,*)amp1,choice(2)
  write(25,fmt="(f8.3,1x,a7,'    # PL-1 amplitude')")amp1,choice(2)

  print*,' Input the exponent 1 < n < 3 and fix/float and Max allowed exponent; default is 4'
  read(5,*,IOSTAT=IOstatus)exp1,choice(3),expmax
  if (IOstatus .ne. 0) expmax=4
  write(25,fmt="(f8.3,1x,a7,2x,f6.2,'    # PL-1 index')")exp1,choice(3),expmax

  if (amp1 .eq. -99999.0) then

!   Make a guess at the colored noise amplitude
     if (rms .gt. est_white_noise) then
        sig_color=sqrt(rms**2 - est_white_noise**2)
        xfny=365.25/(2.0*dt_sam)
        xfl=365.25/(t(ic)-t(1))
        if (exp1 .eq. 1.0 ) then
            PPo=(sig_color**2)/(alog(xfny) - alog(xfl))
        else

          eexp=1.0-exp1
          PPo=(sig_color**2)/(eexp*(xfny**eexp - xfl**eexp))
        end if
!   equation from Langbein (2004) #11
           amp1=((2.0*xfny)**(1.0 - exp1/2.0))*PPo /(2.0*((2.0*3.14159)**(-1.0*exp1))  )
           amp1=sqrt(amp1)
!   If using additive noise, revise the value of white noise downward
           if (( ModType .ne. 'n' ) .and. (siginstr .eq. est_white_noise)) then
              xamp=1.0*siginstr*amp1*((dt_sam/365.25)**(exp1/4))
              if (siginstr**2 .gt. xamp ) siginstr=sqrt(siginstr**2 -xamp)
           end if
     else
        amp1=0.5
     end if
     print*,' Using default amplitude of power law ',amp1 
     print*,' Using default white noise ',siginstr
        

  end if

  print*,' Input the time constant alpha in rad/yr and fix/float'
  read(5,*)alpha,choice(4)
  write(25,fmt="(f8.3,1x,a7,'    # PL-1 G-M term, rad/yr')")alpha,choice(4)
  print*,' '
  print*,' Input the parameters for band-passed filtered noise'
  print*,' Input low and high freq stop band in c/yr'
  read(5,*)f1a,f2a
  if (f1a .lt. f2a) then
     flow=f1a
     fhigh=f2a
   else
     flow=f2a
     fhigh=f1a
  end if
  write(25,fmt="(2f15.5, '    #  bandpass filter elements')")flow,fhigh
  print*,' low frequency stop band is ',flow,' c/yr'
  print*,' high frequency stop band is ',fhigh,' c/yr'
  print*,' number of poles between 1 and 4'
  read(5,*)npole
  write(25,fmt="(i5,'    # number of poles')")npole  
  print*,' Input the amplitude and fix/float'
  read(5,*)amp_bp,choice(5)
  write(25,fmt="(f8.3,1x,a7,'   #  BP amplitude')") amp_bp,choice(5)
  print*,' Input the exponent of second Power Law function fix/float'
  read(5,*)exp2,choice(7)
  write(25,fmt="(f8.3,1x,a7,'    # PL-2 index')")exp2,choice(7)
  print*,' Input the amplitude of second PL function fix/float'
  read(5,*)amp2,choice(6)
  write(25,fmt="(f8.3,1x,a7,'    # PL-2 amplitude')")amp2,choice(6)
  print*,' '
  print*,' Sometimes, it may be necessary to add white noise to data so'
  print*, ' that a better estimate of long period PSD parameters can be made.'
  print*,'This is especially true for data that is predominantly power noise'
  print*,'  but is very rare....'
  print*,' Enter value of white noise to be added (nominal it should be 0)'
  read(5,*)wh_add
  write(25,fmt="(f8.3,'   #  additive white noise')")wh_add

!!  done with guessing

!! Start the calculation clock
  call cpu_time(Tstart)
!  print*,'Tstart',Tstart
!  add white noise to data if wh_add > 0
  if (wh_add .gt. 0.0) then
    call srand(iseed)
    do  i=1,ic
      
       d(i)=d(i)+wh_add*sngl(ZBQLNOR(0.0d+0,1.0d+0))

    end do
  end if
  
!
! zero-out covariance or filter

  if (ModType .eq. 'n' ) then
      allocate(covarpl1(int(itlen*1.00),int(itlen*1.00)))
      allocate(covarpl2(int(itlen*1.00),int(itlen*1.00)))
      allocate(covarbp(int(itlen*1.00),int(itlen*1.00)))
      do  i=1,ic
      do  j=1,i
        covarpl1(i,j)=0.
        covarpl1(j,i)=0.
        covarpl2(i,j)=0.
        covarpl2(j,i)=0.
        covarbp(i,j)=0.
        covarbp(j,i)=0.
     end do
     end do

  else
      maxk=int(alog((itlen*1.00))/alog(2.0))+1
      maxk=2**maxk
      allocate(filtpl1(maxk))
      allocate(filtpl2(maxk))
      allocate(filtbp(maxk))
      do  i=1,int(itlen*1.00)
        filtpl1(i)=0.0
        filtpl2(i)=0.0
        filtbp(i)=0.0
      end do   
  end if
  
!  If the exponentials are to be "fixed" (choice()), calculate filters and store
!  matrix once and store
!
  ipl_flag_1=0
  if ((choice(3) .eq. 'fix') .and. (choice(4) .eq. 'fix'))  then
     print*,' Calculating power law covariance for first set'
     if (ModType .eq. 'n' ) then


        call pow_law_cov(t_year,ic,0.0,1.0,exp1,alpha,covarpl1,size(t_year),size(d),dt_sam/365.25d+0)

     else

        call frac_diff(filtpl1,exp1,alpha,sngl(dt_sam/365.25d+0),size(filtpl1))
     end if
     ipl_flag_1=1
  end if

  ipl_flag_2=0
  if ((choice(7) .eq. 'fix'))  then
     print*,' Calculating power law covariance for second set'
     if (ModType .eq. 'n' ) then

        call pow_law_cov(t_year,ic,0.0,1.0,exp2,alpha2,covarpl2,size(t_year),size(d),dt_sam/365.25d+0)
                    
     else
        alpha2=0.
!        print*,exp2,0.0,sngl(dt_sam/365.25d+0), size(filtpl2),dt_sam
        call frac_diff(filtpl2,exp2,0.0,sngl(dt_sam/365.25d+0),size(filtpl2))
                     
     end if
    ipl_flag_2=1
  end if

    ibp_flag=1
!      if (choice(5) .eq. 'float' ) ibp_flag=0
!  call this once and forget about it!
  if (amp_bp .ne. 0) then
    
    if ((ModType .eq. 'n' ) .or. (ModType .eq. 'q')) then
        print*,' calling band_pass_cov'
!        print*,'ic,max_time,max_data',ic,size(filtbp)
        call band_pass_cov(t_year,ic,dt_sam/365.25d+0,flow,fhigh,npole,size(t_year),size(d),1.0,covarbp)
!                                          (corran,jmax,t_small,flow,fhigh,npole,mmax,md,sig5,covar)
                        
    else
        print*,' call band_pass_filt'

        call band_pass_filt(dt_sam/365.25d+0,flow,fhigh,npole,maxk,filtbp,1.0)
!       print*,(ran4(i),i=1,10)

     end if


  end if

!! allocate variables for different model type
! kmax=int((t_year(jmax)-t_year(1))/t_small + 0.5)+1

  kmax=int((t_year(ic)-t_year(1))/(dt_sam/365.25) + 0.5)+1
  if (Modtype .eq. 'f') call CreateModTypeF
  if (Modtype .eq. 'c') call CreateModTypeC
  if ((Modtype .eq. 'q') .or. (ModType .eq. 'n')) call CreateModTypeQ

!  figure-out how many loops to run the ameoba routine
!
  iln=0
  do  i=1,7
    if (choice(i) .eq. 'float') iln=iln+1
  end do
  if (iln .eq. 1) loop=1
  if (iln .eq. 2) loop=2

  if (iln .gt. 2) loop=iln-1
  if (iln .eq. 0) loop=0
  print*,' Number of loops for downhill simplex', loop
  print*,' ModType ', ModType,  '  calling mle'

!
!  loop through the nedler-meade simplex routine
!
  print*,' '
  print*,' list of trial covariance parameters'
  write(6,fmt="('    white noise     PL_1 amp     PL_1 exp      GM freq','       BP amp' &
      ,'      PL_2 amp      PL_2 exp  determinant  chi^2        MLE ','          cpu ')")

!
!  Big Loop to iterate and find optimal noise model -- only execute when loop > 0 --
!    otherwise, jump close to the end and evaluate funmin and output results without searching
!   for optimal noise
 
  iFlagReLoop=0   
  loopOrig=loop   
!  Jump to here if, after dithering the so-called optimal solution, a better solution was found
      
39130 continue
  if (iFlagReLoop .eq. 1 ) loop=1

    if (loop .ne. 0 ) then
      do  iloop=1,loop    !!!  big loop
!  Initial guess
        iswt=0

        if (iloop .eq. 1) then
       
          call mle(fmle,siginstr,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt)

          fmax=fmle
          do  nx=1,8
             amby(nx)=-fmle
             ambp(nx,1)=siginstr
             ambp(nx,2)=amp1
             ambp(nx,3)=exp1
             ambp(nx,4)=alpha
             ambp(nx,5)=amp_bp
             ambp(nx,6)=amp2
             ambp(nx,7)=exp2
             ifloat(nx)=0
          end do  !! nx
        end if   !! iloop=1



        nx=1
        if (choice(1) .ne. 'fix') then
!  dither the white noise component
!    default dither
          dsig=siginstr*vsig(iloop)
!   compare amplitude of white noise with high frequency of power law noise
!   if power-law high freq is greater than that from white noise; modify dither
          fs=365.25/dt_sam
          psd_db=5.0*exp1*alog10(fs)-7.399*exp1-10*alog10(fs)+2.076
          psd_pl1=(amp1**2)*10.0**(psd_db/10.)
          psd_hf_pl1=psd_pl1/((2.0*fs)**exp1)

          psd_db=5.0*exp2*alog10(fs)-7.399*exp2-10*alog10(fs)+2.076
          psd_pl2=(amp2**2)*10.0**(psd_db/10.)
          psd_hf_pl2=psd_pl2/((2.0*fs)**exp2)
          psd_hf=psd_hf_pl1
          if (psd_hf_pl2 .gt. psd_hf_pl1) psd_hf=psd_hf_pl2
          psd_wh=siginstr**2/(2.0*fs)
          sig_min=sqrt(psd_hf*(2.0*fs))
          if (psd_hf .gt. psd_wh) then
            print*,' Power law noise exceeds that of white noise'
            print*,'  modify white noise from ',siginstr,' to',sig_min
          end if
          if ((psd_hf .gt. psd_wh) .and. (vsig(iloop) .gt. 1)) dsig=vsig(iloop)*sig_min
        
          if ((psd_hf .gt. psd_wh) .and. (vsig(iloop) .lt. 1)) dsig=sig_min/vsig(iloop)
          call mle(fmle,dsig,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt)
 
           ifloat(nx)=1
           nx=nx+1
           amby(nx)=-fmle
           ambp(nx,1)=dsig
           ambp(nx,2)=amp1
           ambp(nx,3)=exp1
           ambp(nx,4)=alpha
           ambp(nx,5)=amp_bp
           ambp(nx,6)=amp2
           ambp(nx,7)=exp2

        end if
        if (choice(2) .ne. 'fix') then
!  dither the power law amplitude
          damp1=amp1*vamp1(iloop)
          call mle(fmle,siginstr,damp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt)
          ifloat(nx)=2
          nx=nx+1
          amby(nx)=-fmle
          ambp(nx,1)=siginstr
          ambp(nx,2)=damp1
          ambp(nx,3)=exp1
          ambp(nx,4)=alpha
          ambp(nx,5)=amp_bp
          ambp(nx,6)=amp2
          ambp(nx,7)=exp2         
        end if
        if (choice(3) .ne. 'fix') then
!  dither the first power law exponent
        dexp1=exp1+vexp1(iloop)
        call mle(fmle,siginstr,amp1,dexp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt)
         ifloat(nx)=3
         nx=nx+1
         amby(nx)=-fmle
         ambp(nx,1)=siginstr
         ambp(nx,2)=amp1
         ambp(nx,3)=dexp1
         ambp(nx,4)=alpha
         ambp(nx,5)=amp_bp
         ambp(nx,6)=amp2
         ambp(nx,7)=exp2
      end if  
      if (choice(4) .ne. 'fix') then
!  dither the GM freq
        dalpha=alpha*valpha(iloop)
!  if 2*pi/alpha is less than longest period in data, modify
!
        alpha_2pi=2.0*3.14159*alpha
        flong=2.0/(t_year(ic)-t_year(1))
        if (alpha_2pi .lt. flong) then
           if (valpha(iloop) .gt. 1) dalpha=flong*valpha(iloop)
           if (valpha(iloop) .lt. 1) dalpha=flong/valpha(iloop)
        end if
  
        call mle(fmle,siginstr,amp1,exp1,dalpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt)
         ifloat(nx)=4
         nx=nx+1
         amby(nx)=-fmle
         ambp(nx,1)=siginstr
         ambp(nx,2)=amp1
         ambp(nx,3)=exp1
         ambp(nx,4)=dalpha
         ambp(nx,5)=amp_bp
         ambp(nx,6)=amp2
         ambp(nx,7)=exp2

      end if        
      if (choice(5) .ne. 'fix') then
! dither the Band pass filter amplitude
        damp_bp=amp_bp*vamp_bp(iloop)
!  compare amplitude of BP filtered noise with that from power law noise;
!    change dither if PL noise exceed BP at 0.5*(flow+fhigh)
        freq=0.5*(flow+fhigh)
        fs=365.25/dt_sam
        psd_db=5.0*exp1*alog10(fs)-7.399*exp1-10*alog10(fs)+2.076
        psd_pl1=(amp1**2)*10.0**(psd_db/10.)
        psd_pl1=psd_pl1/(freq**exp1)
        psd_db=5.0*exp2*alog10(fs)-7.399*exp2-10*alog10(fs)+2.076
        psd_pl2=(amp2**2)*10.0**(psd_db/10.)
        psd_pl2=psd_pl2/(freq**exp2)
        psd_pl=psd_pl1
        if (psd_pl2 .gt. psd_pl1) psd_pl=psd_pl2
          
!  evaluate power level of band pass filtered noise
        h=(1.0**2)*((freq/flow)**2) /((1.0+((freq/flow)**2))*(1.0+((freq/fhigh)**2)))
        h=h**npole

        powBP=h*(amp_bp**2)
        if (psd_pl .gt. powBP) then
            amp_mod=sqrt(psd_pl/h)
            print*,' Modify BP amplitude for dithering'
            print*,' at Freq=',freq
            print*,' flow and fhigh ',flow,fhigh
            print*,'   Power law amplitude is ', psd_pl
            print*,'   BP filter amplitude is ', powBP
            print*,'  Change amplitude from',amp_bp,' to', amp_mod
            print*,'   ignore change in amplitude'
            amp_mod=amp_bp
            if (vamp_bp(iloop) .gt. 1) damp_bp=amp_mod*vamp_bp(iloop)
            if (vamp_bp(iloop) .lt. 1) damp_bp=amp_mod/vamp_bp(iloop)
        end if
         call mle(fmle,siginstr,amp1,exp1,alpha,damp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt) 
         ifloat(nx)=5
         nx=nx+1
         amby(nx)=-fmle
         ambp(nx,1)=siginstr
         ambp(nx,2)=amp1
         ambp(nx,3)=exp1
         ambp(nx,4)=alpha
         ambp(nx,5)=damp_bp
         ambp(nx,6)=amp2
         ambp(nx,7)=exp2
      end if          

      if (choice(6) .ne. 'fix') then
!  dither the second power law amplitude
        damp2=amp2*vamp2(iloop)

        call mle(fmle,siginstr,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,damp2,exp2,ModType,iswt) 
         ifloat(nx)=6

         nx=nx+1
         amby(nx)=-fmle
         ambp(nx,1)=siginstr
         ambp(nx,2)=amp1
         ambp(nx,3)=exp1
         ambp(nx,4)=alpha
         ambp(nx,5)=amp_bp
         ambp(nx,6)=damp2
         ambp(nx,7)=exp2

      end if
      if (choice(7) .ne. 'fix') then
!  dither the second power law exponent
        dexp2=exp2+vexp2(iloop)
        call mle(fmle,siginstr,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,dexp2,ModType,iswt) 

         ifloat(nx)=7
         nx=nx+1
         amby(nx)=-fmle
         ambp(nx,1)=siginstr
         ambp(nx,2)=amp1
         ambp(nx,3)=exp1
         ambp(nx,4)=alpha
         ambp(nx,5)=amp_bp
         ambp(nx,6)=amp2
         ambp(nx,7)=dexp2

      end if

      nfloat=nx-1          
      print*,' Initial solutions for Amoeba'
      do  i=1,nfloat+1
         write(6,fmt='(f12.3,2x,7f7.2)')amby(i),(ambp(i,j),j=1,7)
      end do

!  For initial nelder/mead, us a loose tolerance and progressively
!   tighten tolerance
      FTOL=1.0e-06*(2.5*float(loop) - 2.5*float(iloop-1))
      if ( iFlagReLoop .eq. 1) then
        print*, iFlagReLoop, FTOL,loopOrig
        FTOL=1.0e-06*(2.5*float(loopOrig+1) - 2.5*float(loopOrig))
        iFlagReLoop=0
      end if
      FTOL=FTOL/5.0
      print*,' Tolerance of mle ',FTOL
      ITER=0
      call NedlerMead(ambp,amby,8,7,nfloat,FTOL,ITER,npole,flow,fhigh,ModType,ifloat,iswt)

!      call NedlerMead(ambp,amby,8,7,nfloat,FTOL,ITER,
 !    &   npole,flow,fhigh,dt_sam,max_data,max_mod,
!     &    t_year,iswt,modType,A,d,res,ifloat)


!  pick the smallest solution
      small=1.0e+10

      do  k=1,8
        if (amby(k) .lt. small) then
          small=amby(k)
  
          fmax=amby(k)
          siginstr=ambp(k,1)
          amp1=ambp(k,2)
          exp1=ambp(k,3)
          alpha=ambp(k,4)
          amp_bp=ambp(k,5)
          amp2=ambp(k,6)
          exp2=ambp(k,7)
        end if
      end do
      print*,' '
      print*,' Best fitting solutions'
      print*,' MLE= ',-fmax
      print*,' white noise= ',siginstr
      print*,' Bandpass filter amplitude= ',amp_bp
      print*,' power law noise 1'
      print*,'    amplitude= ',amp1
      print*,'    exponent= ',exp1
      print*,'    G-M freq= ',alpha
      print*,' power law noise 2'
      print*,'    amplitude= ',amp2
      print*,'    exponent= ',exp2
      print*,' '


      do  k=1,8
       amby(k)=fmax
       ambp(k,1)=siginstr
       ambp(k,2)=amp1
       ambp(k,3)=exp1
       ambp(k,4)=alpha
       ambp(k,5)=amp_bp
       ambp(k,6)=amp2
       ambp(k,7)=exp2
      end do




      end do  !  index on iloop
    end if     ! nloop .ne. 0
      


! estimate errors and covariance of spectral parameters

      FTOL=FTOL*abs(fmax)*5.0
      dith(1)=0.01
      dith(2)=0.01
      dith(3)=0.05
      dith(4)=0.05
      dith(5)=0.05
      dith(6)=0.05
      dith(7)=0.1
      xn(1)=siginstr
      xn(2)=amp1
      xn(3)=exp1
      xn(4)=alpha
      xn(5)=amp_bp
      xn(6)=amp2
      xn(7)=exp2

3913  continue

! rerun best fit solution but estimate standard error on model parameters
      iswt=0
      call mle(fmle,siginstr,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt) 

      print*,' Start the covariance calculations for noise model'
      do 3911 i=1,7
        print*,' best estimate',xn(i),'  dither',dith(i)
3911  continue

!  do diagonal terms

  nchoice=0
  do  i=1,7
    ierr=0
   if (choice (i) .eq. 'float') then
        nchoice=nchoice+1
3912    continue
        do k=1,7
          xx(k)=xn(k)
        end do
          xx(i)=xn(i)-dith(i)*xn(i)

          if (ierr .eq. 10 ) go to 39155
!  set max number of loops to get dither! Sometimes, I don't get convergence


        ierr=ierr+1
        call mle(fv2,xx(1),xx(2),xx(3),xx(4),xx(5),npole,flow,fhigh,xx(6),xx(7),ModType,iswt) 
!        call funmin(fv2,xx(1),xx(2),xx(3),xx(4),xx(5),
!     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
!     &    iswt,ModType,A,d,res)


          if (fv2+fmax .gt. 2.0*FTOL) then
             fmax=-fv2
             print*,' Found better, optimal model; re-run NelderMead'
             if (i .eq. 1) siginstr=xx(i)
             if (i .eq. 2) amp1=xx(i)
             if (i .eq. 3) exp1=xx(i)
             if (i .eq. 4) alpha=xx(i)
             if (i .eq. 5) amp_bp=xx(i)
             if (i .eq. 6) amp2=xx(i)
             if (i .eq. 7) exp2=xx(i)
             iFlagReLoop=1
             xn(i)=xx(i)
             go to 39130
!   jump out of do loop since we found a better model!
          end if
!   Modify dither
          if (( abs(fv2 + fmax) .lt. FTOL) .and. (dith(i) .lt. 0.01)) then
             dith(i)=1.9*dith(i) 
             print*,' Dither changed to ', dith(i)
             go to 3912
          end if
          if ((abs(fv2 + fmax) .gt. 5.0*FTOL) .and. (dith(i) .gt. 0.001)) then
            dith(i)=0.6*dith(i)
             print*,' Dither changed to ', dith(i)
            go to 3912
          end if

         
39155     continue

          xx(i)=xn(i)+dith(i)*xn(i)
        call mle(fv1,xx(1),xx(2),xx(3),xx(4),xx(5),npole,flow,fhigh,xx(6),xx(7),ModType,iswt) 
!        call funmin(fv1,xx(1),xx(2),xx(3),xx(4),xx(5),
!     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
!     &    iswt,ModType,A,d,res) 
          if (fv1+fmax .gt. 2.0*FTOL) then

             fmax=-fv1
             print*,' Found better, optimal model; re-run NelderMead'
             if (i .eq. 1) siginstr=xx(i)
             if (i .eq. 2) amp1=xx(i)
             if (i .eq. 3) exp1=xx(i)
             if (i .eq. 4) alpha=xx(i)
             if (i .eq. 5) amp_bp=xx(i)
             if (i .eq. 6) amp2=xx(i)
             if (i .eq. 7) exp2=xx(i)
             iFlagReLoop=1
             xn(i)=xx(i)
             go to 39130
!   jump out of do loop since we found a better model!
          end if


           acov(nchoice,nchoice)=  &
                abs( (fv1+fv2+2*fmax)/( (dith(i)*xn(i))**2 ) )
   


      end if

  end do   
  
!  the off-diagonal terms

  ii=0
  do  i=1,7
    if (choice (i) .eq. 'float') then
      ii=ii+1
      jj=0
      do j=1,i
      if (choice (j) .eq. 'float') then
        jj=jj+1
        do  k=1,5
          xx(k)=xn(k)
        end do
        if (i .ne. j)  then
     
!  i not equal to j....got do 4 of these calculations rather than 2!
            xx(i)=xn(i)+dith(i)*xn(i)
            xx(j)=xn(j)+dith(j)*xn(j)
           call mle(fv1,xx(1),xx(2),xx(3),xx(4),xx(5),npole,flow,fhigh,xx(6),xx(7),ModType,iswt) 

            xx(i)=xn(i)+dith(i)*xn(i)
            xx(j)=xn(j)-dith(j)*xn(j)
            call mle(fv2,xx(1),xx(2),xx(3),xx(4),xx(5),npole,flow,fhigh,xx(6),xx(7),ModType,iswt) 
                   
            xx(i)=xn(i)-dith(i)*xn(i)
            xx(j)=xn(j)+dith(j)*xn(j)
            call mle(fv3,xx(1),xx(2),xx(3),xx(4),xx(5),npole,flow,fhigh,xx(6),xx(7),ModType,iswt) 
           
            xx(i)=xn(i)-dith(i)*xn(i)
            xx(j)=xn(j)-dith(j)*xn(j)
            call mle(fv4,xx(1),xx(2),xx(3),xx(4),xx(5),npole,flow,fhigh,xx(6),xx(7),ModType,iswt) 

!          acov(ii,jj)=abs(fv1+fv4-fv2-fv3)/   sqrt((4.*dith(i)*xn(i)*dith(j)*xn(j)))
          acov(ii,jj)=abs(fv1+fv4-fv2-fv3)/   ((4.*dith(i)*xn(i)*dith(j)*xn(j)))
          acov(jj,ii)=acov(ii,jj)
  
        end if
      end if
      end do
      end if
  end do   !! end off diagonal
!   take inverse of acov
  print*,' '
  ifloat=0
  j=0
  do i=1,7
    if (choice(i) .eq. 'float') then
      j=j+1
      ifloat(j)=i
     
    end if
  end do
  jmax=j
  print*,' Inverse covariance matrix'
  write(6,fmt="(7x,7a10)")(nameNoise(ifloat(i)),i=1,jmax)
  do i=1,jmax
      write(6,fmt="(a7,7f10.1)")nameNoise(ifloat(i)),(acov(i,j),j=1,jmax)
  end do
!

  call ssyev('V','L',nchoice,acov,7,dith,wz,23,ierr)
  print*,' The eigenvalues for inverting covariance matrix',(dith(i),i=1,nchoice)

 
  do  i=1,nchoice
      do  j=1,nchoice
        cov(i,j)=0.0
        do  k=nchoice,1,-1
          cov(i,j)=cov(i,j)+acov(i,k)*acov(j,k)/dith(k)
        end do
    end do
  end do
 print*,' '
  print*,' the covariance matrix'
write(6,fmt="(7x,7(2x,a9))")(nameNoise(i),i=1,7)
!  nzz=0
  ii=0      
  do  i=1,7
      if (choice(i) .eq. 'float') then
        ii=ii+1
!        nzz=nzz+1
        dith(i)=sqrt(abs(acov(ii,ii)))
        jj=0
        do  j=1,i
        if (choice(j) .eq. 'float') then
          jj=jj+1
          acov(ii,jj)=cov(ii,jj)/sqrt(cov(ii,ii)*cov(jj,jj))
          acov(jj,ii)=acov(ii,jj)
          else
          acov(i,j)=0.
        end if
        end do
      else
        dith(i)=0.0
      end if
!       if (choice(i) .eq. 'float') & 
 !        write(6,fmt="(a7,7e11.3)")nameNoise(ifloat(i)),(cov(i,j),j=1,jmax)
 !     print*,i,' cov= ',(cov(i,j),j=1,7)
  end do
  do i=1,jmax
    write(6,fmt="(a7,7e11.3)")nameNoise(ifloat(i)),(cov(i,j),j=1,jmax)
  end do
  print*,' '
  print*,' Cross correlation matrix'
  write(6,fmt="(7x,7a10)")(nameNoise(ifloat(i)),i=1,jmax)
  do  i=1,jmax
!     print*,i,' cross correlation ',(acov(i,j),j=1,7)
     write(6,fmt="(a7,7f10.3)")nameNoise(ifloat(i)),(acov(i,j),j=1,jmax)
  end do 
  if (loop .eq. 0 ) then
! rerun best fit solution but estimate standard error on model parameters
        iswt=0
        call mle(fmle,siginstr,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt) 
 

        fmax=fmle
  end if

! print-out optimal LS fit parameters

!  search 'max.out' for MLE
  print*,'MLE best is ', fmax
  rewind (12)
  rewind (13)
  irow=1
  fmin=9999.e+30
  fmax=-fmax
7245  continue
     read(12,7112,end=7246)fxx
7112    format(91x,2x,13x,13x,f13.3)
  dif=abs(fxx-fmax)
     if (dif .le. fmin) isave=irow
     if (dif .le. fmin) fmin=dif
     irow=irow+1

      go to 7245
7246  continue
  close (12)
  open(12,position="APPEND",file="max.out")

  print*,' '
  print*,' '

  print*,' row number for optimal solution is',isave 
 do i=1,isave-1
        read(13,*)(x(j),ex(j),j=1,nmod)
  end do
  read(13,*)(x(j),ex(j),j=1,nmod)
  print*,' nmod=',nmod
  write(6,fmt="(3x,120(1x,f10.3,' +/- ',f10.3))")(x(j),ex(j),j=1,nmod)

  close(13)
  open(13,position="APPEND",file="model.out")
!  use best MLE and recalculate the model
       iswt=1
!  calculate error bars on model  (x and e) e being the error



!  Apr 2021 --- if alpha .ne. 0  set modtype = 'c'
 if ((alpha .ne. 0 ) .and. ( modtype .eq. 'f') ) then
    Modtype='c'


   call DestructModTypeF
   call CreateModTypeC
!   remove entries in A, d, t_year, and t that were inserted with 0 earlier for 'f'
   k=0
   do i=1,ic
     if (A(i,1) .ne. 0) then
       k=k+1
       do j=1,nmod
         A(k,j)=A(i,j)
       end do
       d(k)=d(i)
       t_year(k)=t_year(i)
       t(k)=t(i)
     end if
   end do
  ic=k
  nmiss=0
 end if

 call mle(fmle,siginstr,amp1,exp1,alpha,amp_bp,npole,flow,fhigh,amp2,exp2,ModType,iswt) 
 if (ModType .eq. 'c') call DestructModTypeC
!  count rows
  irow=0
  rewind(13)
8430   continue
    read(13,*,end=8431)
    
    irow=irow+1
    go to 8430
8431   continue
  print*,' number of rows is', irow
  rewind (13)
  do  i=1,irow-1
       read (13,*)enull
  end do 
  ix=0
  do i=1,n_exp
      if (exp_choice(i) .eq. "float") ix=ix+2
  end do
  read(13,*)(x(j),ex(j),j=1,nmod+ix)
  print*,(x(j),ex(j),j=1,nmod+ix)


!  modify rates due to normalization
  if (nrate .ne. 0) then
    x(nbt+1)=x(nbt+1)/rate_norm
     ex(nbt+1)=ex(nbt+1)/rate_norm
  end if
  if (n_rat_chg .ne. 0) then
    do  kkk=1,n_rat_chg
        x(nbt+nrate+kkk)=x(nbt+nrate+kkk)/rat_chng_norm(kkk)
        ex(nbt+nrate+kkk)=ex(nbt+nrate+kkk)/rat_chng_norm(kkk)
    end do
  end if
  if (n_file_press .ne. 0) then
     nc=nbt+nrate+n_rat_chg+nper*2+noff+nexp_fix
         
     do kkk=1,n_file_press
         x(nc+kkk)=x(nc+kkk)/aux_norm(kkk)
         ex(nc+kkk)=ex(nc+kkk)/aux_norm(kkk)
      end do
  end if

  write(14,fmt="(a20,3x,120(1x,f10.3,' +/- ',f10.3))")filename(1:20),(x(j),ex(j),j=1,nmod),(x(j+nmod),ex(j+nmod),j=1,ix)


!  undo normalization
  if (nrate .ne. 0) then
    x(nbt+1)=x(nbt+1)*rate_norm
    ex(nbt+1)=ex(nbt+1)*rate_norm
  end if
  if (n_rat_chg .ne. 0) then
    do  kkk=1,n_rat_chg
        x(nbt+nrate+kkk)=x(nbt+nrate+kkk)*rat_chng_norm(kkk)
        ex(nbt+nrate+kkk)=ex(nbt+nrate+kkk)*rat_chng_norm(kkk)
    end do
  end if
  if (n_file_press .ne. 0) then
    nc=nbt+nrate+n_rat_chg+nper*2+noff+nexp_fix

    do  kkk=1,n_file_press
         x(nc+kkk)=x(nc+kkk)*aux_norm(kkk)
         ex(nc+kkk)=ex(nc+kkk)*aux_norm(kkk)
    end do
  end if
  
!     nmod=nmod_orig
!  Output the residuals
!


!   JL -- Jan 29, 2016 ---  Noted that A_orig matrix not updated
!    when there is an exponential time constant to be estimated and dec .ne. 0 -- Legacy mode
!   Needed to change ic (which is in common) to add the exponential to A_orig
 !     icTmp=ic
 !     ic=ic_orig
 !     A=A_orig
 !     t=time_orig
      
 !     ic=ic-nmiss
 
      call modify_A
      
 !     call modify_A(A_orig,dtime_orig)
 !     ic=icTmp
 !     print*,' Update A_orig matrix'
!  output model parameters using best noise model
       

      open (65,file='resid.out')
      call outResid
!      call outResid(ic_orig,nmod,A_orig,d_orig,x,dtime_orig,t_start,
!     &  netx,max_data,max_mod)
     
  do i=1,nmod+(n_exp-nexp_fix)

        ex(i)=ex(i)*Anor(i)
  end do      
  print*,' '
  call output
  print*,' '
  print*,' '

  print*," "
  close (65)
  nc=0
  do i=1,7
    ex(i)=0.
    if (choice(i) .eq. 'float') then
        nc=nc+1
        ex(i)=sqrt(abs(cov(nc,nc)))
    end if
  end do
  
  print*,' AIC= ', 2*(nc+nmod+ix) - 2*fmax
  print*,' BIC= ', (nc+nmod+ix)*alog(float(ic_orig)) - 2.0*fmax
  print*,' Best fitting solutions'
  print*,' MLE= ',fmax
  write(6,fmt="('  white noise= ',f12.5, ' +/-', f10.5 )") siginstr,ex(1)
  write(6,fmt="('  Bandpass filter amplitude= ',f12.5, ' +/-', f10.5)") amp_bp,ex(5)
  print*,' power law noise 1'
  write(6,fmt="('    amplitude= ',f12.5, ' +/-', f10.5)")amp1,ex(2)
  write(6,fmt="('    exponent= ',f12.5, ' +/-', f10.5)")exp1,ex(3)
  write(6,fmt="('    G-M freq= ',f12.5, ' +/-', f10.5,' rad/yr')")alpha,ex(4)
  print*,' power law noise 2'
  write(6,fmt="('    amplitude= ',f12.5, ' +/-', f10.5)")amp2,ex(6)
  write(6,fmt="('    exponent= ',f12.5, ' +/-', f10.5)")exp2,ex(7)  
  print*,' '
  icx=ic
! for fast covariance analysis, need to change the actual number of points in max.dat
  if ( ModType .eq. 'f' ) icx=ic-nmiss
  write(15,fmt="(a20,3x,f10.3,3i6,7(3x,f10.3,1x,f10.3),3x,f10.3)   ")filename(1:20),fmax,idec,icx,nmod, &
     siginstr,ex(1),amp1,ex(2),exp1, ex(3),alpha, ex(4), amp_bp, ex(5), amp2, ex(6), exp2, ex(7)
 
  print*,' Residuals in resid.out and/or resid_dec.out'
  print*,'   time, obs-calc, calc, obs'
  print*,' Model covariance and cross correlation in covar.out'
  print*,' Journal of input parameter in estin.jrn'
  print*,' Parameters needed for program adjust_1 in adj.out'
  print*,'   adj.out time format is doy'
  print*,' History of estimating exp/Omori time constants,',' if requested in tauexp.out'
!  compute wall clock time and cpu cycles in second
  call system_clock(count1, count_rate, count_max)
  call cpu_time(timex1)
  write(6,fmt="(' runtime per wall clock',f12.2,' sec')")(dble(count1)-dble(count0))/dble(count_rate)
  print*,' runtime per cpu cycles', timex1-timex0,' seconds'
  print*,' Number of threads used',OMP_get_max_threads()

  if (Modtype .eq. 'f') call DestructModTypeF     
  if (Modtype .eq. 'c') call DestructModTypeC    
  if ((Modtype .eq. 'q') .or. (ModType .eq. 'n')) call DestructModTypeQ    
  call DestructGlobal
  deallocate(A_orig)
  deallocate(time_orig)
  deallocate(d_orig)
  deallocate(t_year_orig)
  if (allocated(covarpl1)) deallocate(covarpl1)
  if (allocated(covarpl2)) deallocate(covarpl2)
  if (allocated(covarbp)) deallocate(covarbp)
  if (allocated(filtpl1)) deallocate(filtpl1)
  if (allocated(filtpl2)) deallocate(filtpl2)
  if (allocated(filtbp)) deallocate(filtbp)
  call random_seed()
  call RANDOM_NUMBER(dist)   !! dist is just a declared variable with no intrinsic meaning
  iseed=int(1000000.*dist)
  rewind (23)
  write(23,*)iseed
!!  Close files
  close(25)
  close(23)
  close(17)
  close(89)
  close(12)
  close(13)
  close(14)
  close(15)
end program est_noise8
