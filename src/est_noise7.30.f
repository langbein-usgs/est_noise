C  Revise est_noise6.50 and est_noise7.10 to combine both Max Likelihood estimator (MLE)
c   to determine noise models to time-series data and Simultaneously, fit a simple temporal
c  function to the time-series.
c
c  Noise model for est_noise6.50 adds several independent noise sources together as
c
c    N(t)^2 = (h_1 * e_1)^2 + (h_2 * e_2)^2 + (h_3 * e_3)^2 ...
c
c   where * is convolution, h_i are impluse responce of a i_th noise model and e_i are
c     independent white noise sources
c
c  An alternate noise model consists of simply adding the impluse responses and then convolving with
c   a single white noise source;
c
c   N(t)=  [ h_1 + h_2 + h_3 +...] * e(t)
c
c   which is used in est_noise7.10 and its variants
c
c   The additive form as the advantage of rapid compution of the inverse of the data
c   covariance, and for missing data, using the method of Bos et al (2012) to speed up
c   the inversion.
c
c   Currently, noise models are
c   Power-law noise as generated using Hosking (1981)
c   Gauss Markove noise discussed by Langbein (2004)
c   white noise
c   Band Pass filtered noise, Langbein (2004)
c   Program could be expanded to ARIMA models using Hosking (1981)
c
c  Temporal dependence includes
c   Nominal value including one or more time-series
c   rate
c   changes in rate over a specified time
c   offsets at prescribed times
c   sinusoids (in phase and quadrature terms)
c   exponential trends
c   Omori law (log(1+t/tau)
c   User prescibed external functions (Pressure data for strainmeter adjustment)
c
c  Program is configured such that standard input used historically for est_noise6.50 will
c   work as default for this new version.
c
c  Program includes option to create noisy data sampled at time of real data for simulations
c
c  Program includes additional time-stamps other than those used befoer;
c   year day-of-year
c   year, month, day (numbers)
c   GMT  (Generic Mappping Tools time stamps)
c   MJD (Modified Julian Day)
c
c  Command to compile on Mac OS 10.6 with Intel compiler and MKL libraries
c  ifort   -Wl,-stack_size,0x40000000 -O3      -L/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal -lmkl_lapack -lmkl_core -lmkl_intel_lp64 -lmkl_intel_ilp64 -lmkl_intel_thread -o est_noise7.22 est_noise7.22.f ~/proglib/sublib/time.f funmin7.22.f invfft.f NedlerMeadSimplex.f  modify.f NedlerMeadSimplexEXP7.22.f
c
c  Using gfortran
c  gfortran -O3 -o est_noise7.21 est_noise7.21.f time.f  funmin7.21.f invfft.f  -framework Accelerate NedlerMeadSimplex.f  modify.f NedlerMeadSimplexEXP7.21.f
c
c  Modify Apr 2021 to correctly compute covariance matriculates when encountering GGM noise model
c
      character*3 net, netx, netp
      character*1 ans,ModType
      character*7 choice(10),exp_choice(10),exp_type(10)
      character*132 filename,string
c  Time is double precision!
      double precision  t_start, t_stop
      double precision t(12687),dtime(12687),
     & dtime_orig(12687),txx,tpress,dt_sam,ZBQLNOR,fjul,
     & t_orig(12687),
     & t_year(12687),t_year_orig(12687),timex,dts
      double precision rat_chg1(82),rat_chg2(82),rat_chng_norm(82),
     &  per(82),off(82),texp(82),aux_norm(82)
      double precision covar(12687,12687), covinv(12687,12687)
      double precision A(12687,82),A_orig(12687,82),
     & fmle,fv1,fv2,fv3,fv4,fmax
      double precision filtpl1(32768),filtpl2(32768),filtbp(32768),
     & covarpl1(12687,12687),covarpl2(12687,12687),
     & covarbp(12687,12687)
      dimension d(12687),d_orig(12687),x(82),e(82),res(12687)
c   ,enull(82)
      dimension bexp(82)
      integer irowmiss(12687)
      dimension  ambp(8,7),amby(8),ifloat(8), dith(7), xn(7),xx(7)
      dimension vsig(8),vamp1(8),vexp1(8),vamp_bp(8),valpha(8),vamp2(8),
     & vexp2(8),acov(7,7),cov(7,7),wz(23)
      double precision ran4(32768),chi2,amby
c   double precision
      common /ModFit1/dtime,t_year,covinv,covar,texp
      common /ModFit1a/t_start,t_stop
c  single precision
      common /ModFit2/x,e,bexp
      common /ModFit2a/expmax
c integer
      common /ModFit3/max_data,max_mod,ic,nmod,n_exp,nmiss,irowmiss,
     & max_time,ipl_flag_1,ipl_flag_2,ibp_flag,nmiss_max
c  character
      common /ModFit4/ exp_choice,exp_type
c  Intermediate covariance and filter functions; double precision
      common /CovFlt1/filtpl1,filtpl2,filtbp,covarpl1,covarpl2,covarbp


      data (vsig(i),i=1,8)/0.7,1.2,1.2,1.2,0.7,1.2,1.2,0.7/
      data (vamp1(i),i=1,8)/1.95,1.95,0.3,0.3,0.3,0.3,1.95,1.95/
      data (vexp1(i),i=1,8)/0.4,-0.4,-0.4,-0.4,-0.4,-0.4,0.4,-0.4/
      data (valpha(i),i=1,8)/4.5,0.35,4.5,4.5,4.5,4.5,0.35,4.5/
      data (vamp_bp(i),i=1,8)/0.55,1.79,0.55,1.79,0.55,1.79,0.55,0.55/
      data (vamp2(i),i=1,8)/2.31,0.26,2.31,0.26,0.26,0.26,2.25,0.26/
      data (vexp2(i),i=1,8)/-0.6,-0.6,0.6,0.6,-0.6,-0.6,0.6,0.6/
      max_time=32768
      max_parm=82
      max_data=12687
      nmiss_max=4103
      max_mod=82
       print*,' Program estimates the power spectral density (PSD)'
       print*,'  of data.  PSD functions can be a combination of:'
       print*,'  1) white noise'
       print*,'   2) Simple power law noise P/f^n'
       print*,'  3) a Gauss Markov version P/(fa^n + f^2)'
       print*,'      where fa=alpha/2*pi'
       print*,'  4) A second power law'
       print*,'  5) Band-pass filtered noise'
       print*,' '
       print*,' While estimating the PSD function, program also',
     &     ' estimates'
       print*,'  various parameters that describe the time series',
     &     ' including'
       print*,'  1) DC term'
       print*,'  2) rate'
       print*,'  3) sinusoidal amplitudes of specified frequencies'
       print*,'  4) rate changes'
       print*,'  5) offset'
       print*,'  6) simple exponential and log(t) function'
       print*,'  7) user supplied function or data'
       print*,' '
       print*,' Program can handle data in various formats'
       print*,' Version 7.30  --- May 2021'
       print*,' '
c  next statement only for MKL/Intel linked libraries
c      print*,' num of threads',mkl_get_max_threads()
      open(23,file='seed.dat')
      open(12,file='max.out')
      open(13,file='model.out')
      open(14,file='model.dat')
      open(15,file='max.dat')
      open(25,file='estin.jrn')
      open(26,file='prob1.out')
      open(78,file='covar.out')
      open(89,file='tauexp.out')
      read(23,*)iseed
      print*,' A journal file consisting of typed input'
      print*,'  is created in estin.jrn'

      print*,'Input the data type for  processing'

      print*,' otr=data'
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

      if ((net .ne. 'mjd') 
     & .and.  (net .ne. 'otr') .and.
     &  (net .ne. 'otd') .and. (net .ne. 'gmt') 
     &   .and. (net .ne. 'otx')) then
        print*,' code for network not in list'
        print*,'Must be either otr otx otd mjd or gmt'
        stop
      end if
      write(25,*)net
c
c  Read Inputs that specifies the time-dependent model
c
      nmod=0
      print*,' input number of time series'
      print*,'   program will estimate only one set of PSD functions'
      read(5,*)nbt
      write(25,301)nbt
301   format(i4, '    # number of time series')
      nmod=nmod+nbt
      print*,' input the period of interest'

      call GetTime1(net,t_start,t_stop)
      print*,' Day numbers from 1960 ',t_start,t_stop
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
      write(25,302)n_rat_chg
302   format(i4,'    #  Number of rate changes')
      nmod=nmod+n_rat_chg
      if (nmod .gt. max_parm) then
         print*,' Number of parmeters exceed limit of', max_parm
         stop
      end if
      if (n_rat_chg .ne. 0) then
        do 9 i=1,n_rat_chg
          call GetTime2(net,i,rat_chg1(i), rat_chg2(i))
          if (rat_chg1(i) .ge. t_stop) then
             print*,' rate change beyond end of time series'
             stop
          end if
          rat_chng_norm(i)=(rat_chg2(i)-rat_chg1(i))/365.25
          print*,' Interval of rate change', i
          print*,rat_chg1(i), rat_chg2(i)
          print*,' Normalization of rate change is:', rat_chng_norm(i)
9       continue
      end if
      print*,' '
      print*,' Input the number of periodicities in the data'
      read(5,*)nper
      write(25,303)nper
303   format(i4,'    # Number of periodicities')
      nmod=nmod+nper*2
      if (nmod .gt. max_parm) then
         print*,' Number of parmeters exceed limit of', max_parm
         stop
      end if
      if (nper .ne. 0) then
         do 10 i=1,nper
         print*,' Input the period (in days) of ', i,' period'
         read(5,*)per(i)
         write(25,*)per(i)
10       continue
      end if
      print*,' '
      print*,' Offset:  specify first day after offset'
      print*,' Input the number of offsets'
      read(5,*)noff
      write(25,304)noff
304   format(i4,'    # Number of offsets')
      nmod=nmod+noff
      if (noff .gt. max_parm) then
         print*,' Number of parmeters exceed limit of', max_parm
         stop
      end if
      if (noff .ne. 0) then
        do 12 i=1,noff
          call GetTime3(net,i,off(i))
          print*,' Offset at ',off(i)
          if (off(i) .gt. t_stop) then
            print*,' Offset beyond end of time series'
            stop
         end if
12    continue
      end if
c  
c  
      print*,' '
      print*,' Input the number of exponentials in time series'
      read(5,*) n_exp
      write(25,305)n_exp
305   format(i4,'    # Number of exponentials/log func.')
      nmod=nmod+n_exp
      nexp_fix=0
      if (n_exp .ne. 0) then
        do 13 i=1,n_exp
          call GetTime4(net,i,texp(i))
          print*,' Exponential starting at ',texp(i)
          if (texp(i) .ge. t_stop) then
             print*,' exponential beyond end of time series'
             stop
          end if     
          print*,' Input the time constant in years and fix/float'
          read(5,*)bexp(i),exp_choice(i)
          write(25,*) bexp(i),exp_choice(i)
          print*,' Input the type of function; e for 1-exp(-t/tau)'
c          print*,'      or                     o for log10(tau + t)'
          print*,'      or                     m for log10(1.0 + t/tau)'
          read(5,*)exp_type(i)
          write(25,*)exp_type(i)
13      continue
        nexp_fix=0
        do  14 k=1,n_exp
        if (exp_choice(k) .eq. 'fix    ') nexp_fix=nexp_fix+1
14      continue
        print*,' The number of exponential time constants that',
     & ' are fixed:', nexp_fix

      end if
c
c  input the data and construc the A matrix
c
      ic=0
      do 50 i=1,nbt
        print*,' Input the format style of baseline data',
     &   ' (otr, otx, otd, mjd, or gmt)'
        read(5,*)netx
        write(25,*)netx
        print*,' Input name of file for baseline number ',i
        read(5,1230)filename
        write(25,*)filename
1230    format(a132)
        open (21,file=filename)
21      continue
        read(21,1230,end=49)string
c        print*,'netx ',netx
c        print*,ic,nmod
c        print*,string
        call GetData(netx,string,timex,dist,err)
c        print*," timex,dist",timex,dist
c  is time within range of t_start and t_stop?
        if ((timex .lt. t_start) .or. (timex .gt. t_stop)) go to 21
c  acceptable data at this point
        ic=ic+1
        dtime(ic)=timex-t_start
        t(ic)=timex-t_start
        dtime_orig(ic)=dtime(ic)
        d(ic)=dist
c
c  build A matrix
c
c  treat different baselines like offsets
        nc=0
        A(ic,i)=1
        nc=nbt+nc
c  secular rate
        if (nrate .eq. 1) A(ic,nc+1)=(t(ic)/365.25)/rate_norm
         nc=nc+nrate
c  rate changes
        if (n_rat_chg .gt. 0) then
          do 25 k=1,n_rat_chg
            A(ic,nc+k)=0.0
            if ((timex .ge. rat_chg1(k)) .and. 
     &       (timex .lt. rat_chg2(k))) A(ic,nc+k)=
     &       ((timex-rat_chg1(k))/365.25)/rat_chng_norm(k)
            if (timex .ge. rat_chg2(k)) A(ic,nc+k)=
     &       ((rat_chg2(k)-rat_chg1(k))/365.25)/rat_chng_norm(k)
25        continue
          end if
          nc=nc+n_rat_chg
c periodicities
        if (nper .gt. 0) then
          do 26 k=1,nper

          txx=2.0*3.141592653589793d0*t(ic)/per(k)
          A(ic,nc+2*k-1)=cos(txx)
          A(ic,nc+2*k-0)=sin(txx)
26        continue
        end if
        nc=nc+2*nper
c  offsets in data
        if (noff .gt. 0) then
          do 27 k=1,noff
          A(ic,nc+k)=0.0
          if (timex .ge. off(k)) A(ic,nc+k)=1
27        continue
        end if
        nc=nc+noff
c
c  Exponentials
c
       n_exp_fix=0
       if (n_exp .ne. 0 ) then
       kz=0
         do 28 k=1,n_exp
         if ( exp_choice(k) .eq. 'fix') then
         kz=kz+1
         A(ic,nc+kz)=0.0
         if (timex .gt. texp(k)) then
           if (exp_type(k) .eq. "e")
     &      A(ic,nc+kz)=1.0-exp(-(timex-texp(k))/(bexp(k)*365.25))
           if (exp_type(k) .eq. "o")
     &     A(ic,nc+kz)=alog10(abs(bexp(k)) +  
     &     (sngl(timex-texp(k))/365.25))
           if (exp_type(k) .eq. "m")
     &     A(ic,nc+kz)=alog10(1.0 +
     &       (sngl(timex-texp(k))/365.25)/abs(bexp(k)) )
         end if
         end if
28       continue
       nexp_fix=kz
       end if
       nc=nc+nexp_fix
        write(26,101)ic,t(ic),d(ic), (A(ic,j),j=1,nc)
101       format(i5,f15.5,f10.2,40(3x,e13.4))
        go to 21
49      continue
        close(21)
50    continue
c
c  input the pressure data and add to A matrix
c    Ugly code!!
c
c
      print*,' column of A matrix prior to input pressure data', nc
      print*,' Number of files of Auxillary data (pressure)'
      read(5,*)n_file_press
      write(25,306)n_file_press
306   format(i4,'    # Number of auxillary time series')
      if (n_file_press .gt. 0) then
        open(88,file='dist_press.dat')
        do 40 i=1,n_file_press

          print*,' Input the format style of baseline data',
     &   ' (otr, otx, otd)'
          read(5,*)netp
          write(25,*)netp
          print*,' Input name of file for aux data number ',i
          read(5,1230)filename
          write(25,*)filename
          open (22,file=filename)
          print*,' Input the lag in days'
          read(5,*)t_lag
          write(25,*)t_lag
          press_min=999999.
          press_max=-999999.
          sum_press=0.0
          fmiss=-9999.0e+20
c  put missing data symbol in pressure column of A matrix
          do 441 ik=1,ic
441       A(ik,nc+i)=fmiss
44        continue
c  read pressure data
          iflag=0
         read(22,1230,end=48)string
         call GetData(netp,string,timex,fpress,err)
         tpress=timex-t_lag

         small=1.0e-04
           tpress=tpress-t_start

         if ((tpress .lt. dtime(1)-0.001) .or. 
     &          (tpress .gt. dtime(ic)+0.001)) then
           go to 44
           else
c      time of pressure in range of data

c    scan through data until tpress=t(izz)
           izz=1
41         continue
c           diff_time=t(izz)-tpress
           diff_time=dtime(izz)-tpress
           if (abs(diff_time) .le. 0.0007) then
c  THIS is a Match!
c              write(88,*)' Match achieved   izz=',izz,' fpress=',fpress

              A(izz,nc+i)=fpress
              if ((fpress) .gt. press_max)
     &          press_max=fpress
              if ((fpress) .lt. press_min)
     &           press_min=fpress
              sum_press=sum_press+fpress
c      flag is set to go read next pressure measurement
              iflag=1
              
            end if
            if (diff_time .lt. -0.0007) then
c  no patch, time of observation less than time of pressure measurement
c        increment time of observation
              izz=izz+1
              iflag=0
            end if
            if ( (diff_time .gt. 0.0007)) then
 
c  No Match; time of observation exceeds time of pressure data


c  go read the next pressure data
              iflag=1
              
            end if
c  determine whether to read the next pressure measurement (iflag=1)
c     or increment the time of the data (iflag=0)
            if (iflag .eq. 1) go to 44
            if (iflag .eq. 0) go to 41

         end if


48    continue
c  all pressure data has been read (end of record)
      close (22)

c  scan A matrix for missing pressure data
       ixx=1
443    continue
          if (A(ixx,nc+i) .eq. fmiss) then
               write(88,*) ' no pressure ixx=',ixx,' t(ixx)=',t(ixx),
     &         'ic-1',ic-1

c  toss-out measurement and shuffle a matrix
               do 42 ik=ixx,ic-1
               d(ik)=d(ik+1)
               t(ik)=t(ik+1)
               if (ik .eq. ixx) then
c                 write(88,*)'ixx=',ixx,' t(ixx)=',t(ixx),a(ik,nc+i),
c     &      'ik+1=',ik+1,a(ik+1,nc+i)
               end if
                do 43 jk=1,nc+i
                A(ik,jk)=A(ik+1,jk)
43              continue
42             continue
              ic=ic-1
              ixx=ixx-1
          end if
          ixx=ixx+1
        if (ixx .le. ic) go to 443
       


      aux_norm(i)=(press_max-press_min)
c      aux_norm(i)=1
      ave_press=sum_press/float(ic)
      print*,' Renormalize auxilary data set ',i,' by',aux_norm(i)
      print*,' max press is ', press_max,'   Min press is: ',
     &  press_min,'   Average press is: ', ave_press

      do 47 ik=1,ic
c       write(89,*)t(ik),d(ik),a(ik,nc+i)
       write(88,88890)ik,t(ik),d(ik),a(ik,nc+i),
     &   (a(ik,nc+i)-ave_press)/aux_norm(i)
88890   format(i9,f20.9,3f20.3)
      a(ik,nc+i)=(a(ik,nc+i)-ave_press)/aux_norm(i)
47    continue

      print*,' List of distance vs pressure data ',
     &  'in dist_press.dat'



40    continue
      close (88)
      close (89)


      end if
c  end input of auxillary data
c
      nmod=nbt+nrate+n_rat_chg+2*nper+noff+n_file_press+nexp_fix
c      print*,nbt,nrate,n_rat_chg,2*nper,noff,n_file_press,nexp_fix
      print*,' Number of data read is ',ic
      print*,' Number of model parameters is ',nmod
      

c123   format(f15.5,f8.2,30(1x,f8.5))
      if (nmod .gt. max_mod) then
        print*,' Number of model parameters', nmod,
     &  ' exceeds dimension ',max_mod
        stop
      end if
      if (ic .gt. max_data) then
        print*,' Number of  data', ic,
     &  ' exceeds dimension ',max_data
        stop
      end if
c
c  put data and A in cronological order
c
      call chron(max_data,max_mod,ic,nmod,t,dtime,d,A)
      open(7,file='prob2.out')
      do 52 i=1,ic
52    write(7,101)i,t(i),d(i),(A(i,j),j=1,nmod)
      close (7)
      print*, ' '
      print*,'  prob1.out; time, data, and A matrix ',
     &  ' aux data excluded'

      print*,'  prob2.out has stuff in cronological and aux. data '
c
c  get statistics on sampling interval
c
      dt_ave=(t(ic)-t(1))/float(ic-1)
      dt_min=999999.
      dt_max=-99999.
      dt_sd=0
      do 53 i=2,ic
        if (t(i)-t(i-1) .gt. dt_max) dt_max=t(i)-t(i-1)
        if (t(i)-t(i-1) .lt. dt_min) dt_min=t(i)-t(i-1)
        dt_sd=dt_sd+ (t(i)-t(i-1) -dt_ave)**2
53    continue
      dt_sd=sqrt(dt_sd/float(ic-2))
      print*, ' '
      print*,' Average Sampling interval:  ',dt_ave, ' days'
      print*,' Standard deviation of Ave. interval: ',dt_sd,' days'
      print*,' Shortest Sampling interval: ',dt_min, ' days'
      print*,' Longest Sampling Interval:  ',dt_max, ' days'
      print*,' '
      print*,' Input the minimum sampling interval in days to use'
      read(5,*)dt_sam
      write(25,307)dt_sam
307   format(f20.15,  '    # sampling inteval in days')
      if (dt_min .lt. 0.1*dt_ave) then
      print*," The following may be near simultaneous observations"
      print*,"   and on might need be eliminated"
      do 5301 i=2,ic
        if ( t(i)-t(i-1) .le. 0.1*dt_ave) then
          write(6,5302)i-1,i,t(i-1),d(i-1),t(i),d(i),t(i)-t(i-1) 
        end if
5301  continue
      end if
5302  format("i-1 and i",2i5," t(i-1),d(i-1)",f18.9,f9.2,
     &  " t(i),d(i)",f18.9,f9.2,"  dt=",f12.9)
c
c examine measurements times and if there are redundent measurements
c   within period of less than the minimum sampling interval, average
c   those measurements
c
      npts=int((t(ic)-t(1))/dt_sam)+1
      t_len=t(ic)-t(1)
      if (npts .gt. max_time) then
        print*,' Length, in time, of time series is ',npts,
     &   ' exceeds dimension of ',max_time
        stop
      end if
      print*,' Number of points in time series for analysis is',npts
      print*,' '
      print*,' Input the type Noise model construction'
      print*,' n or q = noise components are summed as squares'
      print*,'   as in est_noise6.x (Langbein, 2004)'
      print*,' a = noise component are simple sums'
      print*,'   automatically switches between cholesky decomposition'
      print*,'    or computing inverse using combo of inverse',
     &    '    convoluton and Bos et al (2012) for missing data'
      print*,' c = noise component simple sums'
      print*,'     forces inversion using cholesky decomposition'
      print*,' f = noise component simple sums'
      print*,'    forces combo of inverse convolution and',
     &  '     and Bos et al (2012) for missing data'
      print*,' '
      read(5,*)ModType
      write(25,*)ModType
      if ( ModType .eq. 'q') ModType='n'
      print*,' Do you want to substitute real data with random',
     &  ' numbers? y/n'
      read(5,*)ans
      write(25,*)ans
      if (ans .eq. 'y') then
c
c  substitute noise for real data; noisy data has same sampling interval as real data
c   NOTE -- This only works correctly if a single baseline is being examined! (nbt=1)
c
        open(77,file='noise.dat')
        call genNoise(ModType,t,d,ic,max_data,iseed,
     &     dt_sam,max_time)
        do 58 i=1,ic
          itime=int(t(i)+t_start)
          dec_time=(t(i)+t_start)-int(t(i)+t_start)
          call inv_jul_time(itime,nyr,jul)
          fjul=float(jul)+dec_time
          write(77,*)nyr,fjul,d(i),t(i)/365.25
 
58      continue
       
        print*,' random noise data written to noise.dat'
        close(77)
      end if
c
c save original data and A matrix
c  estimate high frequency standard error in data by taking 1st differences
c

      ic_orig=ic
      sum_diff=0.
      n_diff=0
      do 601 i=1,ic
        t_year(i)=(t(i)-t(1))/365.25d+00
        t_year_orig(i)=t_year(i)
        d_orig(i)=d(i)
        t_orig(i)=t(i)
        if (i .ne. 1) then
        if (t(i) - t(i-1) .le. dt_sam) then
           sum_diff=sum_diff+(d(i)-d(i-1))**2
           n_diff=n_diff+1
        end if
        end if
        do 602 j=1,nmod
602     A_orig(i,j)=A(i,j)
601   continue
      est_white_noise=sqrt(sum_diff/float(n_diff))/sqrt(2.0)
      print*,' Estimate of white noise component of data is:',
     &  est_white_noise
c
c  decimate the data!  ---  IMPLIMENTED  only for ModType=n 
c
      print*,' methods to decimate data'
      print*,' 0 = no decimation'
      print*,' 1 = keep 2, skip 1, keep 2, skip 1'
      print*,' 2 = keep 2, skip 1, keep 1, skip 2, keep 2, skip 1...'
      print*,' 3 = keep 2, skip 1, keep 1, skip 2, keep 1, skip 3,',
     &  ' keep 2 skip 1...'
      print*,' 4 or more...make it option 3'
      print*,'input choice, BUT is noise model is simple sums,',
     &   ' then ignore'
      read(5,*)idec
      if (ModType .ne. 'n') idec=0
      if (idec .ge. 4) idec=3
      write(25,308)idec
308   format(i4,'     # Decimation type')
      if (idec .ne. 0) then

        call decimate(max_data,max_mod,ic,nmod,A,t,dtime,t_year,d,idec)
        open (77,file='decimate.dat')
        do 59 i=1,ic
          itime=int(t(i)+t_start)
          dec_time=(t(i)+t_start)-int(t(i)+t_start)
          call inv_jul_time(itime,nyr,jul)
          fjul=float(jul)+dec_time
          write(77,*)nyr,fjul,d(i),t(i)/365.25
59      continue
        close(77)
        print*,' Decimated data in decimate.dat'        
      end if
c
c  Estimate the data completeness
c
      nptsThr=int((dtime_orig(ic)-dtime_orig(1))/dt_sam +0.005)
      nmissEst=nptsThr-ic+1
      perNmissEst=100.*(float(nmissEst)/float(nptsThr))
      print*,nptsThr,ic,nmissEst
      print*,' Estimated number of missing data is ', nmissEst
      print*,' Estimated percenta of data gaps is', perNmissEst
        if ( ModType .eq. 'a' ) then
          ModType='f'
          if (perNmissEst .gt. 25.0) ModType='c'
        print*,' ModType changed from a-automatic, to ',ModType
        end if


c
c fit data using design matrix and ASSUMPTION of white noise
c
      do 61 i=1,ic
      covar(i,i)=1.0d+0
      covinv(i,i)=1.0d+0
61    continue
      nopt=1


      call model_fit(nopt,'n',A,d,res,chi2)

      rms=sqrt(chi2/float(ic-nmod))
      print*,' RMS fit using white noise model is',
     &  rms
      det=float(ic-nmiss)*alog(rms)
      f_mle=-det-0.5*chi2/(rms**2)
     &   -0.5*float(ic)*alog(2.0*3.14159265)
      print*,' Log MLE for white noise error model is ',f_mle
      print*, ' '
c  modify white noise version of error bars by scaling by rms
      do i=1,nmod
        e(i)=e(i)*rms
      end do

c
c  output model parameters assuming white noise
c

      call output(x,e,max_mod,net,nbt,nrate,n_rat_chg,
     & nper,noff,rat_chg1,rat_chg2,off,per,sea,rate_norm,
     &  n_file_press,aux_norm,n_exp,texp,bexp,exp_choice,exp_type,
     &  rat_chng_norm)


c
c  For ModType  simple sums (ModType .ne. n)
c  Create a revised design matrix (and data) that fills in data-gaps.
c     force rows of revised design matrix to be filled with zeros
c

      nmiss=0
      if (ModType .eq. 'f') then
      do 201 i=1,ic
        if (i .gt. 1) then
          idiff=int(((dtime_orig(i)-dtime_orig(i-1))/dt_sam)+0.005)
          if (idiff .gt. 1 ) then
c            print*,"i,idiff,dtime(i),dtime(i-1)",i,
c     &       idiff,dtime_orig(i),dtime_orig(i-1)

            do 202 k=1,idiff-1
              nmiss=nmiss+1
              ix=int(((dtime_orig(i)-dtime_orig(1))/dt_sam)+0.005)+1
c              irowmiss(nmiss)=ix-(idiff-1)+k-1+nmiss0
              irowmiss(nmiss)=ix-(idiff-1)+k-1
c              print*,"i,ix,nmiss,irowmiss",i,ix,nmiss,irowmiss(nmiss)
              d(i+nmiss-1)=0.
              dtime(i+nmiss-1)=dtime_orig(i-1)+float(k)*dt_sam
              t_year(i+nmiss-1)=t_year_orig(i-1)
     &          +float(k)*dt_sam/365.25d+0
              do 204 j=1,nmod
                A(i+nmiss-1,j)=0.
204           continue
202         continue
            d(i+nmiss)=d_orig(i)
            dtime(i+nmiss)=dtime_orig(i)
            t_year(i+nmiss)=t_year_orig(i)
            do 203 j=1,nmod
203         A(i+nmiss,j)=A_orig(i,j)       
          else
            d(i+nmiss)=d_orig(i)
            dtime(i+nmiss)=dtime_orig(i)
            t_year(i+nmiss)=t_year_orig(i)
            do 205 j=1,nmod
205         A(i+nmiss,j)=A_orig(i,j)      
          end if
        end if
201   continue
      print*,' Number of missing data are', nmiss
      ic=ic+nmiss
      print*,' percentage of data gaps', 100*(float(nmiss)/float(ic))

        
      end if



      print*,' sampling interval in yrs ',dt_sam/365.25
      print*,' '
c      if (f_mle .lt. 9999999.) stop
      sig1=rms
      print*,'  Input the initial parameters of the PSD and whether',
     &       '   the item is "fix" or "float"'
      print*,'   if "fix", then the item is not estimated'
      print*,'   if "float", then the item is estimated'
      print*,' '
      print*,' Input the white noise "instrument precision" ',
     &          'and fix/float'
      print*,'  Suggested input is ', est_white_noise,
     &  ' enter -99999.0 to use suggested value'
      read(5,*)siginstr, choice(1)
      write(25,309) siginstr, choice(1)
      if (siginstr .eq. -99999.0 ) siginstr=est_white_noise
309   format(f8.3,1x,a7,'    # white noise')
      print*,' Input the amplitude first Power law function ',
     &   'and fix/float'
      print*,' To use default value, enter -99999.0'
      read(5,*)amp1,choice(2)
      write(25,310)amp1,choice(2)
310   format(f8.3,1x,a7,'    # PL-1 amplitude')

      print*,' Input the exponent 1 < n < 3 and fix/float',
     &  ' and Max allowed exponent; default is 4'
      read(5,*,IOSTAT=IOstatus)exp1,choice(3),expmax
      if (IOstatus .ne. 0) expmax=4
      write(25,311)exp1,choice(3),expmax

      if (amp1 .eq. -99999.0) then

c   Make a guess at the colored noise amplitude
        if (rms .gt. est_white_noise) then
           sig_color=sqrt(rms**2 - est_white_noise**2)
c           print*,rms,est_white_noise,sig_color
           xfny=365.25/(2.0*dt_sam)
           xfl=365.25/t_len
c           print*,xfny,xfl,t_len
           if (exp1 .eq. 1.0 ) then
             PPo=(sig_color**2)/(alog(xfny) - alog(xfl))
c             print*,PPo
           else

             eexp=1.0-exp1
             PPo=(sig_color**2)/(eexp*(xfny**eexp - xfl**eexp))
c             print*,PPo
           end if
c   equation from Langbein (2004) #11
           amp1=((2.0*xfny)**(1.0 - exp1/2.0))*PPo
     &         /(2.0*((2.0*3.14159)**(-1.0*exp1))  )
           amp1=sqrt(amp1)
c   If using additive noise, revise the value of white noise downward
           if (( ModType .ne. 'n' ) .and. 
     &          (siginstr .eq. est_white_noise)) then
               xamp=1.0*siginstr*amp1*((dt_sam/365.25)**(exp1/4))
c               print*,siginstr**2,xamp,amp1,dt_sam
               if (siginstr**2 .gt. xamp ) 
     &           siginstr=sqrt(siginstr**2 -xamp)

           end if
        else
          amp1=0.5
        end if
        print*,' Using default amplitude of power law ',amp1 
        print*,' Using default white noise ',siginstr
        

      end if
311   format(f8.3,1x,a7,2x,f6.2'    # PL-1 index')
      print*,' Input the time constant alpha in c/yr and fix/float'
      read(5,*)alpha,choice(4)
      write(25,312)alpha,choice(4)
312   format(f8.3,1x,a7,'    # PL-1 G-M term, c/yr')
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
      write(25,313)flow,fhigh
313   format(2f15.5, '    #  bandpass filter elements')
      print*,' low frequency stop band is ',flow,' c/yr'
      print*,' high frequency stop band is ',fhigh,' c/yr'
      print*,' number of poles between 1 and 4'
      read(5,*)npole
      write(25,314)npole
314   format(i5,'    # number of poles')
      print*,' Input the amplitude and fix/float'
      read(5,*)amp_bp,choice(5)
      write(25,315) amp_bp,choice(5)
315   format(f8.3,1x,a7,'   #  BP amplitude')
      print*,' Input the exponent of second Power Law function',
     &  ' fix/float'
      read(5,*)exp2,choice(7)
      write(25,316)exp2,choice(7)
316   format(f8.3,1x,a7,'    # PL-2 index')
      print*,' Input the amplitude of second PL function fix/float'
      read(5,*)amp2,choice(6)
      write(25,317)amp2,choice(6)
317   format(f8.3,1x,a7,'    # PL-2 amplitude')
      print*,' '
      print*,' Sometimes, it may be necessary to add white noise to ',
     &'data so'
      print*, ' that a better estimate of long period PSD parameters ',
     &'can be made.'
      print*,'This is especially true for data that is predominantly ',
     & 'power noise'
      print*,' Enter value of white noise to be added (nominal it ',
     &'should be 0)'
      read(5,*)wh_add
      write(25,318)wh_add
318   format(f8.3,'   #  additive white noise')

c  Start the calculation clock
      Tstart=secnds(0.)
c  add white noise to data if wh_add > 0
      if (wh_add .gt. 0.0) then
        call srand(iseed)
        do 8889 i=1,ic
        d(i)=d(i)+wh_add*sngl(ZBQLNOR(0.0d+0,1.0d+0))

8889    continue
      end if
c
c zero-out covariance or filter
c
      if (ModType .eq. 'n' ) then
      do 69 i=1,ic
      do 691 j=1,i
        covarpl1(i,j)=0.
        covarpl1(j,i)=0.
        covarpl2(i,j)=0.
        covarpl2(j,i)=0.
        covarbp(i,j)=0.
        covarbp(j,i)=0.
691   continue
69    continue

      else

      do 692 i=1,max_time
      filtpl1(i)=0.0
      filtpl2(i)=0.0
      filtbp(i)=0.0
692   continue    
      end if
c
c  If the exponentials are to be "fixed" (choice()), calculate filters and store
c  matrix once and store
c
      ipl_flag_1=0
      if ((choice(3) .eq. 'fix') .and. (choice(4) .eq. 'fix'))  then
         print*,' Calculating power law covariance for first set'
         if (ModType .eq. 'n' ) then


           call pow_law_cov(t_year,ic,0.0,1.0,exp1,alpha,
     &       covarpl1,max_time,max_data,dt_sam/365.25d+0)

         else

           call frac_diff(filtpl1,exp1,alpha,sngl(dt_sam/365.25d+0),
     &       max_time,
     &       max_time)
         end if
         ipl_flag_1=1
      end if

      ipl_flag_2=0
      if ((choice(7) .eq. 'fix'))  then
         print*,' Calculating power law covariance for second set'
         if (ModType .eq. 'n' ) then
           call pow_law_cov(t_year,ic,0.0,1.0,exp2,0.0,
     &       covarpl2,max_time,max_data,dt_sam/365.25d+0)

         else
           call frac_diff(filtpl2,exp2,0.0,sngl(dt_sam/365.25d+0),
     &       max_time,
     &        max_time)
         end if
         ipl_flag_2=1
      end if

      ibp_flag=1
c      if (choice(5) .eq. 'float' ) ibp_flag=0
c  call this once and forget about it!
      if (amp_bp .ne. 0) then
        if (ModType .eq. 'n' ) then
           print*,' calling band_pass_cov'
           print*,'ic,max_time,max_data',ic,max_time,max_data
           call band_pass_cov(t_year,ic,dt_sam/365.25d+0,flow,
     &       fhigh,npole,max_time,max_data,1.0,covarbp)
        else
          print*,' call band_pass_filt'
          call band_pass_filt(dt_sam/365.25d+0,flow,fhigh,npole,
     &      max_time,max_time/2,ran4,1.0)
c       print*,(ran4(i),i=1,10)
          do 693 i=1,max_time/2

c          if (abs(ran4(i)) .lt. 1.0e-08) ran4(i)=1.0e-08
C        write(71,*)float(i)/365,ran4(i)
            filtbp(i)=ran4(i)
693       continue
        end if


      end if
c
c  figure-out how many loops to run the ameoba routine
c
      iln=0
      do 3999 i=1,7
      if (choice(i) .eq. 'float') iln=iln+1
3999  continue
      if (iln .eq. 1) loop=1
      if (iln .eq. 2) loop=2
c      if (iln .gt. 2) loop=iln
      if (iln .gt. 2) loop=iln-1
      if (iln .eq. 0) loop=0
      print*,' Number of loops for downhill simplex', loop
      print*,' ModType ', ModType,  '  calling funmin'
c
c  Dicy exercise of interpolating the data and A matrix to fill in missing data form ModType=f
c
c      if ( ModType .eq. 'f') then
c        print*,' interpolate'
c        do 555 i=1,nmiss
c          print*,i,irowmiss(i)
c          d0=d(irowmiss(i)-1)
c          if (irowmiss(i+1) -irowmiss(i) .gt. 1) then
c            d1=d(irowmiss(i)+1)
c            d(irowmiss)=0.5*(d1 + d0)
c            print*,d(irowmiss(i)),d1,d0
c          end if
           

c555     continue
c      end if

c
c  loop through the ameoba routine
c
      print*,' '
      print*,' list of trial covariance parameters'
      write(6,6666)
6666  format('    white noise     PL_1 amp     PL_1 exp      GM freq',
     &'       BP amp',
     &'      PL_2 amp      PL_2 exp  determinant  chi^2        MLE ',
     &'          cpu ')






c
c  Big Loop to iterate and find optimal noise model -- only execute when loop > 0 --
c    otherwise, jump close to the end and evaluate funmin and output results without searching
c    for optimal noise
c
      iFlagReLoop=0   
      loopOrig=loop   
c  Jump to here if, after dithering the so-called optimal solution, a better solution was found
      
39130 continue
      if (iFlagReLoop .eq. 1 ) loop=1

      if (loop .ne. 0 ) then
      do 801 iloop=1,loop
c  Initial guess
      iswt=0

      if (iloop .eq. 1) then

      call funmin(fmle,siginstr,amp1,exp1,alpha,amp_bp,
     &   npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &  iswt,ModType,A,d,res)

c       print*,'iloop ',iloop
c       stop
c       if (iloop .lt. 100) stop

       fmax=fmle
       do 802 nx=1,8
       amby(nx)=-fmle
       ambp(nx,1)=siginstr
       ambp(nx,2)=amp1
       ambp(nx,3)=exp1
       ambp(nx,4)=alpha
       ambp(nx,5)=amp_bp
       ambp(nx,6)=amp2
       ambp(nx,7)=exp2
       ifloat(nx)=0
802    continue
       end if
       nx=1
      if (choice(1) .ne. 'fix') then
c  dither the white noise component
c    default dither
        dsig=siginstr*vsig(iloop)
c   compare amplitude of white noise with high frequency of power law noise
c   if power-law high freq is greater than that from white noise; modify dither
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
          print*,'  modify white noise from ',siginstr,' to',
     &    sig_min
        end if
        if ((psd_hf .gt. psd_wh) .and. (vsig(iloop) .gt. 1))
     &    dsig=vsig(iloop)*sig_min
        
        if ((psd_hf .gt. psd_wh) .and. (vsig(iloop) .lt. 1))
     &    dsig=sig_min/vsig(iloop)
        call funmin(fmle,dsig,amp1,exp1,alpha,amp_bp,
     &   npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)
 
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
c  dither the power law amplitude
        damp1=amp1*vamp1(iloop)
c  compare with amplitude of band-pass filtered noise at the lowest periods
c    modify dither if PL amplitude is less that BP filtered noise
        if (choice(5) .ne. 'fix') then
          freq=2.0/(t_year(ic)-t_year(1))
c  evalutate  power PL noise
          fs=365.25/dt_sam
          psd_db=5.0*exp1*alog10(fs)-7.399*exp1-10*alog10(fs)+2.076
          psd_pl1=(amp1**2)*10.0**(psd_db/10.)
          psd_pl=psd_pl1/(freq**exp1)
c  evaluate power level of band pass filtered noise
          h=(1.0**2)*((freq/flow)**2)
     &      /((1.0+((freq/flow)**2))*(1.0+((freq/fhigh)**2)))
          if (npole .eq. 2) h=h*h
          if (npole .eq. 3) h=h*h*h
          if (npole .eq. 4) h=h*h*h*h
          if (npole .eq. 5) h=h*h*h*h*h
          if (npole .eq. 6) h=h*h*h*h*h*h
          powBP=h*(amp_bp**2)
          if (powBP .gt. psd_pl) then
            psd_pl=powBP*(freq**exp1)
            amp_mod=psd_pl/(10.**(psd_db/10))
            amp_mod=sqrt(amp_mod)
            if (vamp1(iloop) .gt. 1) damp1=amp_mod*vamp1(iloop)
            if (vamp1(iloop) .lt. 1) damp1=amp_mod/vamp1(iloop)
            print*,'BP Noise at freq=', freq,' exceeds PL noise'
            print*,' PL amp at freq is',psd_pl,' BP amp at freq is',
     &        powBP
            print*,'  Change PL scale from', amp1,' to',amp_mod
          end if
        end if
        call funmin(fmle,siginstr,damp1,exp1,alpha,amp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  

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
c  dither the power law exponent
        dexp1=exp1+vexp1(iloop)

        call funmin(fmle,siginstr,amp1,dexp1,alpha,amp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  
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
c  dither the GM freq
        dalpha=alpha*valpha(iloop)
c  if 2*pi/alpha is less than longest period in data, modify
c
        alpha_2pi=2.0*3.14159*alpha
        flong=2.0/(t_year(ic)-t_year(1))
        if (alpha_2pi .lt. flong) then
           if (valpha(iloop) .gt. 1) dalpha=flong*valpha(iloop)
           if (valpha(iloop) .lt. 1) dalpha=flong/valpha(iloop)
        end if
  
        call funmin(fmle,siginstr,amp1,exp1,dalpha,amp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  
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
c  dither the Band pass filter amplitude
        damp_bp=amp_bp*vamp_bp(iloop)
c  compare amplitude of BP filtered noise with that from power law noise;
c    change dither if PL noise exceed BP at 0.5*(flow+fhigh)
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
          
c  evaluate power level of band pass filtered noise
          h=(1.0**2)*((freq/flow)**2)
     &      /((1.0+((freq/flow)**2))*(1.0+((freq/fhigh)**2)))
          if (npole .eq. 2) h=h*h
          if (npole .eq. 3) h=h*h*h
          if (npole .eq. 4) h=h*h*h*h
          if (npole .eq. 5) h=h*h*h*h*h
          if (npole .eq. 6) h=h*h*h*h*h*h
          powBP=h*(amp_bp**2)
          if (psd_pl .gt. powBP) then
            amp_mod=sqrt(psd_pl/h)
            print*,' Modify BP amplitude for dithering'
            print*,' at Freq=',freq
            print*,'   Power law amplitude is ', psd_pl
            print*,'   BP filter amplitude is ', powBP
            print*,'  Change amplitude from',amp_bp,' to', amp_mod
            if (vamp_bp(iloop) .gt. 1) damp_bp=amp_mod*vamp_bp(iloop)
            if (vamp_bp(iloop) .lt. 1) damp_bp=amp_mod/vamp_bp(iloop)
          end if
    
        call funmin(fmle,siginstr,amp1,exp1,alpha,damp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  
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
c  dither the second power law amplitude
        damp2=amp2*vamp2(iloop)

        call funmin(fmle,siginstr,amp1,exp1,alpha,amp_bp,
     &    npole,flow,fhigh,damp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  

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
c  dither the second power law exponent
        dexp2=exp2+vexp2(iloop)
 
        call funmin(fmle,siginstr,amp1,exp1,alpha,amp_bp,
     &    npole,flow,fhigh,amp2,dexp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)   
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
c      do 888 i=1,nfloat
c      print*," ifloat(i)",i,ifloat(i)
c888   continue
      print*,' Initial solutions for Amoeba'
      do 4030 i=1,nfloat+1
4030  write(*,'(f12.3,2x,7f7.2)'),amby(i),(ambp(i,j),j=1,7)
      FTOL=0.002
      FTOL=0.005
      FTOL=5.0e-06*(1.0 - 0.5*float(iloop-1)/float(loop+1))
c  For initial nelder/mead, us a loose tolerance and progressively
c   tighten tolerance
      FTOL=1.0e-06*(2.5*float(loop) - 2.5*float(iloop-1))
      if ( iFlagReLoop .eq. 1) then
        print*, iFlagReLoop, FTOL,loopOrig
        FTOL=1.0e-06*(2.5*float(loopOrig+1) - 2.5*float(loopOrig))
        iFlagReLoop=0
      end if
      FTOL=FTOL/5.0
      print*,' Tolerance of mle ',FTOL
      ITER=0
c      call AMOEBA(ambp,amby,8,7,nfloat,FTOL,ITER,
       call NedlerMead(ambp,amby,8,7,nfloat,FTOL,ITER,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat)

c  pick the smallest solution
      small=1.0e+10
      do 769 k=1,8
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
769   continue
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


      do 770 k=1,8
       amby(k)=fmax
       ambp(k,1)=siginstr
       ambp(k,2)=amp1
       ambp(k,3)=exp1
       ambp(k,4)=alpha
       ambp(k,5)=amp_bp
       ambp(k,6)=amp2
       ambp(k,7)=exp2
770    continue




801   continue
C  End Big Loop for iterations



c
c estimate errors and covariance of spectral parameters
c
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



c rerun best fit solution but estimate standard error on model parameters
      iswt=0

        call funmin(fmle,siginstr,amp1,exp1,alpha,amp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)   

      print*,' Start the covariance calculations for noise model'
      do 3911 i=1,7
        print*,' best estimate',xn(i),'  dither',dith(i)
3911  continue

c
c  do diagonal terms
c
      nchoice=0
      do 3915 i=1,7
        ierr=0
        if (choice (i) .eq. 'float') then
        nchoice=nchoice+1
3912    continue
        do 3914 k=1,7
3914    xx(k)=xn(k)
          xx(i)=xn(i)-dith(i)*xn(i)
c  comment-out if xx(i) <0
c          if (xx(i) .lt. 0) then
c            xx(i)=0.
c            dith(i)=1.0
c          end if
c
          if (ierr .eq. 10 ) go to 39155
c  set max number of loops to get dither! Sometimes, I don't get convergence
c

        ierr=ierr+1
        call funmin(fv2,xx(1),xx(2),xx(3),xx(4),xx(5),
     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)

c          if (fv2+fmax .gt. 10.0*FTOL) then
c          if (fv2+fmax .gt. 5.0*FTOL) then
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
c   jump out of do loop since we found a better model!
          end if
c   Modify dither
          if (( abs(fv2 + fmax) .lt. FTOL) .and.
     &      (dith(i) .lt. 0.01)) then
             dith(i)=1.9*dith(i) 
             print*,' Dither changed to ', dith(i)
             go to 3912
          end if
          if ((abs(fv2 + fmax) .gt. 5.0*FTOL) .and.
     &     (dith(i) .gt. 0.001)) then
            dith(i)=0.6*dith(i)
             print*,' Dither changed to ', dith(i)
            go to 3912
          end if

         
39155     continue

          xx(i)=xn(i)+dith(i)*xn(i)

        call funmin(fv1,xx(1),xx(2),xx(3),xx(4),xx(5),
     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res) 
          if (fv1+fmax .gt. 2.0*FTOL) then
c          if (fv1+fmax .gt. 5.0*FTOL) then
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
c   jump out of do loop since we found a better model!
          end if
c          if (iflag .eq. 1) i=5
c          if (iflag .eq. 1) go to 3913
           acov(nchoice,nchoice)=
     &     abs( (fv1+fv2-2*fmax)/( (dith(i)*xn(i))**2 ) )

c          print*,' test ',nchoice,nchoice,fv1,fv2,fmax,dith(i)*xn(i),
c     &     acov(nchoice,nchoice),dith(i),xn(i)
c          print*,(fv1+fv2-2*fmax),( (dith(i)*xn(i))**2 )
      end if

3915  continue

c  the off-diagonal terms
c
      ii=0
      do 3920 i=1,7
      if (choice (i) .eq. 'float') then
      ii=ii+1
      jj=0
      do 3919 j=1,i
      if (choice (j) .eq. 'float') then
        jj=jj+1
        do 3902 k=1,5
3902    xx(k)=xn(k)
        if (i .ne. j)  then
c     
c  i not equal to j....got do 4 of these calculations rather than 2!
          xx(i)=xn(i)+dith(i)*xn(i)
          xx(j)=xn(j)+dith(j)*xn(j)

        call funmin(fv1,xx(1),xx(2),xx(3),xx(4),xx(5),
     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)

          xx(i)=xn(i)+dith(i)*xn(i)
          xx(j)=xn(j)-dith(j)*xn(j)
 
        call funmin(fv2,xx(1),xx(2),xx(3),xx(4),xx(5),
     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)

          xx(i)=xn(i)-dith(i)*xn(i)
          xx(j)=xn(j)+dith(j)*xn(j)

        call funmin(fv3,xx(1),xx(2),xx(3),xx(4),xx(5),
     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)

          xx(i)=xn(i)-dith(i)*xn(i)
          xx(j)=xn(j)-dith(j)*xn(j)

        call funmin(fv4,xx(1),xx(2),xx(3),xx(4),xx(5),
     &    npole,flow,fhigh,xx(6),xx(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)
          acov(ii,jj)=abs(fv1+fv4-fv2-fv3)/
     &        (4.*dith(i)*xn(i)*dith(j)*xn(j))
          acov(jj,ii)=acov(ii,jj)
c          print*,ii,jj,
c     &     fv1,fv2,fv3,fv4,dith(i)*x(i),dith(j)*x(j),acov(i,j),
c     &     fv1+fv4-fv2-fv3,(4.*dith(i)*xn(i)*dith(j)*xn(j))
        end if
      end if
3919  continue
      end if
3920  continue

c   take inverse of acov
      print*,' '
      print*,' Inverse covariance matrix'
      do 3930 i=1,7
      print*,i, 'acov= ',(acov(i,j),j=1,7)
3930  continue
c  invert acov  to get covariance matrix

      call ssyev('V','L',nchoice,acov,7,dith,wz,23,ier)
      print*,' The eigenvalues ',(dith(i),i=1,nchoice)

      print*,' ier=', ier
      do 3931 i=1,nchoice
      do 3932 j=1,nchoice
      at=0.
        do 3933 k=nchoice,1,-1
        at=at+acov(i,k)*acov(j,k)/dith(k)
3933    continue
      cov(i,j)=at
3932  continue
3931  continue
      print*,' '
      print*,' the covariance matrix'
      nzz=0
      ii=0      
      do 3935 i=1,7
      if (choice(i) .eq. 'float') then
        ii=ii+1
        nzz=nzz+1
        dith(i)=sqrt(abs(acov(ii,ii)))
        jj=0
        do 3934 j=1,i
        if (choice(j) .eq. 'float') then
          jj=jj+1
          acov(ii,jj)=cov(ii,jj)/sqrt(cov(ii,ii)*cov(jj,jj))
          acov(jj,ii)=acov(ii,jj)
          else
          acov(i,j)=0.
        end if
3934    continue
      else
        dith(i)=0.0
      end if
      print*,i,' cov= ',(cov(i,j),j=1,7)
3935  continue
      print*,' '
      print*,' Cross correlation matrix'
      do 3940 i=1,7
      print*,i,' cross correlation ',(acov(i,j),j=1,7)
3940  continue      

c   end if statement for loop = 0
      end if
      if (loop .eq. 0 ) then
c rerun best fit solution but estimate standard error on model parameters
        iswt=0

        call funmin(fmle,siginstr,amp1,exp1,alpha,amp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)   

        fmax=fmle
      end if
c
c print-out optimal LS fit parameters
c
c  search 'max.out' for MLE
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
      do 7247 i=1,isave
      read(13,*)(x(j),e(j),j=1,nmod)
7247  continue 
      read(13,*)(x(j),e(j),j=1,nmod)
      print*,' nmod=',nmod
      write(6,72488)(x(j),e(j),j=1,nmod)
      print*,(x(j),e(j),j=1,nmod)
72488  format(3x,120(1x,f10.3,' +/- ',f10.3))
c      close (13)
c      open(13,position="APPEND",file="model.out")
c  use best MLE and recalculate the model
       iswt=1
c  calculate error bars on model  (x and e) e being the error

c   Apr 2021 --- if alpha .ne. 0  set modtype = 'c'
       if ((alpha .ne. 0 ) .and. 
     &    ( modtype .eq. 'f') )  Modtype='c'
c     &    ((Modtype .eq. 'f') .or. (Modtype .eq. 'q')) Modtype='c'

        call funmin(fv4,siginstr,amp1,exp1,alpha,amp_bp,
     &    npole,flow,fhigh,amp2,exp2,dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)

         print*
       rewind (13)
c  count rows
       irow=0
8430   continue
         read(13,*,end=8431)
         irow=irow+1
         go to 8430
8431   continue
       print*,' number of rows is', irow
       rewind (13)
       do 8432 i=1,irow-1
       read (13,*)enull
8432   continue
      ix=0
      do 8433 i=1,n_exp
8433  if (exp_choice(i) .eq. "float") ix=ix+2
c      print*," nmod,ix", nmod,ix
      read(13,*)(x(j),e(j),j=1,nmod),(x(j+nmod),e(j+nmod),j=1,ix)
c      read(13,*)(x(j),e(j),j=1,nmod),(x(j+nmod),enull,j=1,ix)
c      write(6,7248)filename(1:20),(x(j),e(j),j=1,nmod),
c     &  (x(j+nmod),e(j+nmod),j=1,ix)
      print*,(x(j),e(j),j=1,nmod),
     &  (x(j+nmod),e(j+nmod),j=1,ix)
c     &  (x(j+nmod),enull,j=1,ix)

c  modify rates due to normalization
      if (nrate .ne. 0) then
        x(nbt+1)=x(nbt+1)/rate_norm
        e(nbt+1)=e(nbt+1)/rate_norm
      end if
      if (n_rat_chg .ne. 0) then
        do 72481 kkk=1,n_rat_chg
        x(nbt+nrate+kkk)=x(nbt+nrate+kkk)/rat_chng_norm(kkk)
        e(nbt+nrate+kkk)=e(nbt+nrate+kkk)/rat_chng_norm(kkk)
72481   continue
      end if
      if (n_file_press .ne. 0) then
         nc=nbt+nrate+n_rat_chg+nper*2+noff+nexp_fit
         
         do 72483 kkk=1,n_file_press
         x(nc+kkk)=x(nc+kkk)/aux_norm(kkk)
         e(nc+kkk)=e(nc+kkk)/aux_norm(kkk)
72483    continue
      end if

      write(14,7248)filename(1:20),(x(j),e(j),j=1,nmod),
     &  (x(j+nmod),e(j+nmod),j=1,ix)

7248  format(a20,3x,120(1x,f10.3,' +/- ',f10.3))
c  undo normalization
      if (nrate .ne. 0) then
        x(nbt+1)=x(nbt+1)*rate_norm
        e(nbt+1)=e(nbt+1)*rate_norm
      end if
      if (n_rat_chg .ne. 0) then
        do 72482 kkk=1,n_rat_chg
        x(nbt+nrate+kkk)=x(nbt+nrate+kkk)*rat_chng_norm(kkk)
        e(nbt+nrate+kkk)=e(nbt+nrate+kkk)*rat_chng_norm(kkk)
72482   continue
      end if
      if (n_file_press .ne. 0) then
         nc=nbt+nrate+n_rat_chg+nper*2+noff+nexp_fit

         do 72484 kkk=1,n_file_press
         x(nc+kkk)=x(nc+kkk)*aux_norm(kkk)
         e(nc+kkk)=e(nc+kkk)*aux_norm(kkk)
72484    continue
      end if
c
c  use optimal model to compute data residuals
c
      nmod_orig=nmod
c
c  NEED TO FIX
c      call modify_A(A,max_data,max_mod,ic,nmod,dtime,
c     &  n_exp,bexp,texp,exp_choice,exp_type,t_start)


c same for undecimated data
      nmod=nmod_orig
c NEED TO FIX
c      call modify_A(A_orig,max_data,max_mod,ic_orig,nmod,dtime_orig,
c     &  n_exp,bexp,texp,exp_choice,exp_type,t_start)

C   JL -- Jan 29, 2016 ---  Noted that A_orig matrix not updated
c    when there is an exponential time constant to be estimated and dec .ne. 0 -- Legacy mode
c   Needed to change ic (which is in common) to add the exponential to A_orig
      icTmp=ic
      ic=ic_orig
      call modify_A(A_orig,dtime_orig,texp)
      ic=icTmp
      print*,' Update A_orig matrix'
c  output model parameters using best noise model
       

      open (65,file='resid.out')
      call outResid(ic_orig,nmod,A_orig,d_orig,x,dtime_orig,t_start,
     &  netx,max_data,max_mod)
      print*," "
      print*,' Residualdata in resid.out'
      print*,' col 1 & 2, time; 3 is residual,',
     &  ' 4 is calculated, 5 is data'
      print*," "
c
      print*,' '
      print*,' '
      call output(x,e,max_mod,net,nbt,nrate,n_rat_chg,
     & nper,noff,rat_chg1,rat_chg2,off,per,sea,rate_norm,
     &  n_file_press,aux_norm,n_exp,texp,bexp,exp_choice,exp_type,
     &  rat_chng_norm)
      print*,' '
      print*,' '

      print*," "
      close (65)
      nc=0
      do 9015 i=1,7
        e(i)=0.
        if (choice(i) .eq. 'float') then
           nc=nc+1
           e(i)=sqrt(abs(cov(nc,nc)))
        end if
9015  continue
c      print*,'nc= ',nc,' nmod= ,',nmod,' ic_orig=',ic_orig,
c     & ' ix= ',ix
      print*,' AIC= ', 2*(nc+nmod+ix) - 2*fmax
      print*,' BIC= ', (nc+nmod+ix)*alog(float(ic_orig)) - 2.0*fmax
      print*,' Best fitting solutions'
      print*,' MLE= ',fmax
      print*,' white noise= ',siginstr,' +/-',e(1)
      print*,' Bandpass filter amplitude= ',amp_bp,
     &    ' +/-',e(5)
      print*,' power law noise 1'
      print*,'    amplitude= ',amp1,' +/-',e(2)
      print*,'    exponent= ',exp1,' +/-',e(3)
      print*,'    G-M freq= ',alpha,' +/-',e(4)
      print*,' power law noise 2'
      print*,'    amplitude= ',amp2,' +/-',e(6)
      print*,'    exponent= ',exp2,' +/-',e(7)
      print*,' '
      icx=ic
c for fast covariance analysis, need to change the actual number of points in max.dat
      if ( ModType .eq. 'f' ) icx=ic-nmiss
      write(15,7879)filename(1:20),fmax,idec,icx,nmod,
     &  siginstr,e(1),
     &  amp1,e(2),    
     &  exp1, e(3),    
     &  alpha, e(4),    
     &  amp_bp, e(5),    
     &  amp2, e(6),    
     &  exp2, e(7),sea
7879  format(a20,3x,f10.3,3i6,7(3x,f10.3,1x,f10.3),3x,f10.3)    
      print*,' Residuals in resid.out and/or resid_dec.out'
      print*,' Model covariance and cross correlation in covar.out'
      print*,' Journal of input parameter in estin.jrn'
      print*,' History of estimating exp/Omori time constants,',
     & ' if requested in tauexp.out'
   
      close(12)
      close(13)
      close(14)
      close(15)
      close(78)
      close(79)
c  Compute execution time
      Tend=secnds(0.)
      TotalTime=Tend-Tstart
      write(6,78791)TotalTime
78791  format(' Total execution time',f20.2," sec")

      iseed=int(1000000.*rand())
      rewind (23)
      write(23,*)iseed
      close(23)

      stop
      end

      subroutine chron(nd,md,ic,nmod,t,dtime,d,a)
c  puts data into chronological order
c  code riped-off from numerical recipes piksrt.f
c
c  nd and md are dimensions of t, d, and a
c  t is time   dtime is same but double precision
c  d is data
c  a is the design matrix
c  ic is number of data
c  nmod number of model parameters
      dimension t(nd),d(nd),a(nd,md),ax(md)
      double precision dtime(nd), dtx,t,a,ax
      do 12 j=2,ic
         tx=t(j)
         dtx=dtime(j)
         dx=d(j)
         do 9 k=1,nmod
9        ax(k)=a(j,k)

         do 11 i=j-1,1,-1
c  if time is a duplicate, add 0.00001 days to value (about 1 second)
         if (t(i) .eq. tx) t(i)=t(i)+0.00001
         if (t(i) .lt. tx) go to 10
           t(i+1)=t(i)
           dtime(i+1)=dtime(i)
           d(i+1)=d(i)
           do 13 k=1,nmod
13         a(i+1,k)=a(i,k)
11       continue
         i=0
10       t(i+1)=tx
         dtime(i+1)=dtx
         d(i+1)=dx
         do 15 k=1,nmod
15       a(i+1,k)=ax(k)
12    continue
      return
      end

      subroutine GetTime1(net,t_start,t_stop)
c
c   read time interval to analyze data
c
c  net is input for the time format
c  t_start and _stop are output: time span, in days since 1960, to analyze data
c
      character*3 net
      character*80 string1,string2
      double precision fjul1, fjul2, t_start, t_stop,day1,day2
      if (net .eq. 'otr') then
        print*,' start and stop times; year day_of_yr year day_of_yr'
        print*,'  "day_of_yr" may be decimal day'
        read(5,*)yr1,fjul1,yr2,fjul2
        write(25,*)yr1,fjul1,yr2,fjul2
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        t_start=dble(iday)+fjul1-int(fjul1)

        call jul2day(int(yr2),int(fjul2),iday)
c days since 1960
        t_stop=dble(iday)+fjul2-int(fjul2)
      end if 
      if (net .eq. 'otx' ) then
        print*,' start and stop times; year mn da  year mn da'
        print*,'   "da" may be a decimal day'
        read(5,*) yr1,mn1,day1,yr2,mn2,day2
        write(25,*) yr1,mn1,day1,yr2,mn2,day2
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        call date2day(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'otd' ) then
        print*,' start and stop times; YearMmnDa  YearMmnDay'
        print*,'   "da" may be a decimal day'
        read(5,*) day1, day2
        write(25,*) day1, day2
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        yr2=float(int(day2/10000.))
        mn2=int((day2-10000*yr2)/100.)
        day2=day2-10000.*yr2- 100.0*mn2
        call date2day(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'mjd' ) then
        print*,' start and stop times; MJD  MJD'
        print*,'   May be decimal days'
        read(5,*) day1, day2
        write(25,*) day1, day2
        t_start=day1-36933.0
        t_stop=day2-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*," start and stop times in GMT time format"
        print*,"   two entries required"
        read(5,*) string1, string2
        write(25,*) string1, string2
        read(string1,100)iyr,imn,ida,ihr,mn,sec
100     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)
        call date2day(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)

        read(string2,100)iyr,imn,ida,ihr,mn,sec
        call date2day(iyr,imn,ida,iday)
        t_stop=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)
      end if
      return
      end
      subroutine GetTime2(net,n,t_start,t_stop)
c
c   read time interval to estimate a RATE change
c
c  net is input for the time format
c  n is an index of the rate change
c  t_start and _stop are output: time span, in days since 1960, to analyze data
c
      character*3 net
      character*80 string1,string2
      double precision fjul1, fjul2, t_start, t_stop,day1,day2
      if (net .eq. 'otr') then
        print*,'rate change interal',n,' year day_of_yr year day_of_yr'
        print*,'  "day_of_yr" may be decimal day'
        read(5,*)yr1,fjul1,yr2,fjul2
        write(25,*)yr1,fjul1,yr2,fjul2
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        t_start=dble(iday)+fjul1-int(fjul1)

        call jul2day(int(yr2),int(fjul2),iday)
c days since 1960
        t_stop=dble(iday)+fjul2-int(fjul2)
      end if 
      if (net .eq. 'otx' ) then
        print*,'  rate change interal',n,' year mn da  year mn da'
        print*,'   "da" may be a decimal day'
        read(5,*) yr1,mn1,day1,yr2,mn2,day2
        write(25,*) yr1,mn1,day1,yr2,mn2,day2
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        call date2day(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'otd' ) then
        print*,' rate change interal',n,' YearMmnDa  YearMmnDay'
        print*,'   "da" may be a decimal day'
        read(5,*) day1, day2
        write(25,*) day1, day2
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        yr2=float(int(day2/10000.))
        mn2=int((day2-10000*yr2)/100.)
        day2=day2-10000.*yr2- 100.0*mn2
        call date2day(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'mjd' ) then
        print*,'  rate change interal',n,' MJD  MJD'
        print*,'   May be decimal days'
        read(5,*) day1, day2
        write(25,*) day1, day2
        t_start=day1-36933.0
        t_stop=day2-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*,'  rate change interal',n,' in GMT time format'
        print*,'   two entries required'
        read(5,*) string1, string2
        write(25,*) string1, string2
        read(string1,100)iyr,imn,ida,ihr,mn,sec
100     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)
        call date2day(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)

        read(string2,100)iyr,imn,ida,ihr,mn,sec
        call date2day(iyr,imn,ida,iday)
        t_stop=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)
      end if
      return
      end
      subroutine GetTime3(net,n,t_start)
c
c   read time to estimate an OFFSETs
c
c  net is input for the time format
c  n is an index of the rate change
c  t_start day since 1960

c
      character*3 net
      character*80 string1
      double precision fjul1, t_start,day1
      if (net .eq. 'otr') then
        print*,' Offset time',n,' year day_of_yr '
        print*,'  "day_of_yr" may be decimal day'
        read(5,*)yr1,fjul1
        write(25,*)yr1,fjul1
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        t_start=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        print*,'  Offset time',n,' year mn da  '
        print*,'   "da" may be a decimal day'
        read(5,*) yr1,mn1,day1
        write(25,*) yr1,mn1,day1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        print*,' Offset Time',n,' YearMmnDa  '
        print*,'   "da" may be a decimal day'
        read(5,*) day1
        write(25,*) day1
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        print*,'  Offset time',n,' MJD '
        print*,'   May be decimal days'
        read(5,*) day1
        write(25,*) day1
        t_start=day1-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*,'  Offset time',n,' in GMT time format'

        read(5,*) string1
        write(25,*) string1
        read(string1,100)iyr,imn,ida,ihr,mn,sec
100     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)
        call date2day(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)
      end if
      return
      end
      subroutine GetTime4(net,n,t_start)
c
c   read time  to estimate an exponential/Omori trend
c
c  net is input for the time format
c  n is an index of the rate change
c  t_start day since 1960
c
      character*3 net
      character*80 string1
      double precision fjul1, t_start,day1
      if (net .eq. 'otr') then
        print*,' Exponential time',n,' year day_of_yr '
        print*,'  "day_of_yr" may be decimal day'
        read(5,*)yr1,fjul1
        write(25,*)yr1,fjul1
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        t_start=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        print*,' Exponential time',n,' year mn da  '
        print*,'   "da" may be a decimal day'
        read(5,*) yr1,mn1,day1
        write(25,*) yr1,mn1,day1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        print*,' Exponential Time',n,' YearMmnDa  '
        print*,'   "da" may be a decimal day'
        read(5,*) day1
        write(25,*) day1
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        print*,' Exponential time',n,' MJD '
        print*,'   May be decimal days'
        read(5,*) day1
        write(25,*) day1
        t_start=day1-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*,' Exponential time',n,' in GMT time format'

        read(5,*) string1
        write(25,*) string1
        read(string1,100)iyr,imn,ida,ihr,mn,sec
100     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)
        call date2day(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)
      end if
      return
      end

      subroutine decimate(max_data,max_parm,ic,nmod,a,t,dtime,
     &   t_year,d,idec)
c idec=1 uses 67% of data
c idec=2 uses 50% of data
c idec=3 uses 40% of data
      dimension t(max_data),d(max_data)
      dimension ax(max_data,max_parm),ttemp(max_data),dtemp(max_data)
      double precision t_year(max_data),tyrtemp(max_data),
     &  dtime(max_data),dtimetmp(max_data),t,ax,ttemp,
     &  a(max_data,max_parm)
      icdec=1
c      print*,idec,d(2),t(2),t_year(2)
      if (idec .eq. 1) then
         do 6501 i=4,ic-4,3
         dtemp(icdec)=d(i-3)
         ttemp(icdec)=t(i-3)
         dtimetmp(icdec)=dtime(i-3)
         tyrtemp(icdec)=t_year(i-3)
           do 6511 j=1,nmod
           ax(icdec,j)=a(i-3,j)
6511       continue
         icdec=icdec+1
         dtemp(icdec)=d(i-2)
         ttemp(icdec)=t(i-2)
         dtimetmp(icdec)=dtime(i-2)
         tyrtemp(icdec)=t_year(i-2)
           do 6521 j=1,nmod
           ax(icdec,j)=a(i-2,j)
6521       continue
         icdec=icdec+1
6501     continue
      end if
      if (idec .eq. 2) then
         do 6502 i=7,ic-7,6
         dtemp(icdec)=d(i-6)
         ttemp(icdec)=t(i-6)
         dtimetmp(icdec)=dtime(i-6)
         tyrtemp(icdec)=t_year(i-6)
           do 6512 j=1,nmod
           ax(icdec,j)=a(i-6,j)
6512       continue
         icdec=icdec+1
         dtemp(icdec)=d(i-5)
         ttemp(icdec)=t(i-5)
         dtimetmp(icdec)=dtime(i-5)
         tyrtemp(icdec)=t_year(i-5)
           do 6522 j=1,nmod
           ax(icdec,j)=a(i-5,j)
6522       continue
         icdec=icdec+1
         dtemp(icdec)=d(i-3)
         ttemp(icdec)=t(i-3)
         dtimetmp(icdec)=dtime(i-3)
         tyrtemp(icdec)=t_year(i-3)
           do 6532 j=1,nmod
           ax(icdec,j)=a(i-3,j)
6532       continue
         icdec=icdec+1
6502     continue
      end if
      if (idec .eq. 3) then
         do 6503 i=11,ic-11,10
         dtemp(icdec)=d(i-10)
         ttemp(icdec)=t(i-10)
         dtimetmp(icdec)=dtime(i-10)
         tyrtemp(icdec)=t_year(i-10)
           do 6513 j=1,nmod
           ax(icdec,j)=a(i-10,j)
6513       continue
         icdec=icdec+1
         dtemp(icdec)=d(i-9)
         ttemp(icdec)=t(i-9)
         dtimetmp(icdec)=dtime(i-9)
         tyrtemp(icdec)=t_year(i-9)
           do 6523 j=1,nmod
           ax(icdec,j)=a(i-9,j)
6523       continue
         icdec=icdec+1
         dtemp(icdec)=d(i-7)
         ttemp(icdec)=t(i-7)
         dtimetmp(icdec)=dtime(i-7)
         tyrtemp(icdec)=t_year(i-7)
           do 6533 j=1,nmod
           ax(icdec,j)=a(i-7,j)
6533       continue
         icdec=icdec+1
         dtemp(icdec)=d(i-4)
         ttemp(icdec)=t(i-4)
         dtimetmp(icdec)=dtime(i-4)
         tyrtemp(icdec)=t_year(i-4)
           do 6543 j=1,nmod
           ax(icdec,j)=a(i-4,j)
6543       continue
         icdec=icdec+1
6503     continue
      end if
         
      if (idec .ne. 0) then
      print*,' Using ',icdec-1,' out of ',ic,' values'
      ic=icdec-1
      do 6599 i=1,ic
      d(i)=dtemp(i)
      t(i)=ttemp(i)
      dtime(i)=dtimetmp(i)
      t_year(i)=tyrtemp(i)
        do 6598 j=1,nmod
6598    a(i,j)=ax(i,j)
6599  continue
      end if
      return
      end
      subroutine model_fit(nopt,ModType,A,d,res,chi2)


c  Subroutine check whether exp_fit or omori_fit has any "floats";
c    if so, used downhill simplex to estimate the optimal time constants
c    if not, just runs calcres
c
c max_data  dimension of data
c max_parm  dimension of num parameters
c ic number of data
c nmod number of parameters
c A  design matrix
c d  data
c t  time
c res  residuals
c covinv  inverse covariance matrix
c covar  covariance matrix
c nopt  outputs the model and its error if set to 1
c  x  the model
c  e  the error
c  nexp  number of exponentials both "fix" and "float"
c  bexp  are the time constants for the exponential (floaters are trial time constrants
c  texp  are the T_o of each of the exponentials
c  exp_choice is either "fix" or "float"
c  t_start  time, in days of first pt in time series
c  t_stop    time, in day of last pt in time series
c  chi2  res^t*covinv*res  chi square
c


      double precision A(max_data,max_mod),dtime(12687),
     & covinv(12687,12687),covar(12687,12687),chi2,sres(11)
      double precision t_start,t_stop,texp(82),sumRes,t_year(12687)
      dimension d(max_data),res(max_data),x(82),e(82),
     & bexp(82),xjunk(20),Anor(82)
      character*1 ModType
      character*7 exp_choice(10),exp_type(10)
      dimension spar(11,10)
      integer irowmiss(12687)
      double precision a1(max_mod,max_data),a2(max_mod,max_mod),
     &  ai(max_mod,max_mod),wz(5000)
      double precision ad2(max_mod,max_mod),wzd(5000),evd(max_mod)
      double precision Finv(32768),Aw(max_data,max_mod),B(12687,12687),
     & din(max_data),dout(max_data),Awt(max_mod,max_data),
     &  a3(max_mod,max_data)

 
c   double precision
      common /ModFit1/dtime,t_year,covinv,covar,texp
      common /ModFit1a/t_start,t_stop
c  single precision
      common /ModFit2/x,e,bexp
      common /ModFit2a/expmax
c integer
      common /ModFit3/max_data,max_mod,ic,nmod,n_exp,nmiss,irowmiss,
     & max_time,ipl_flag_1,ipl_flag_2,ibp_flag,nmiss_max

c  character
      common /ModFit4/ exp_choice,exp_type
c  double precision
      common /Modfit5/ Finv, B

      mdmax=max_data
      maxnmod=max_mod
      nmod_orig=nmod
      nexp=n_exp
      nopt2=0



      if (nexp .ne. 0 ) then


        call modify_A(A,dtime,texp)

      end if

      call calcres(nopt,ModType,sumRes,A,d,res,chi2)

      sres(1)=chi2

c
c  do a second, third, etc solution depending upon the number of "floaters"
c
      nsim=0
      rewind (23)
      read(23,*)iseed
      call srand(iseed)
      nmod=nmod_orig
      iflag=0
      tmin=99999.
      do 111 i=2,ic
        dt=dtime(i)-dtime(i-1)
        if (dt .lt. tmin) tmin=dt
111   continue
      tmax=dtime(ic)-dtime(i)
      kcnt=0
c      print*,' nexp=',nexp
      do 1 k=1,nexp
        if (exp_choice(k) .eq. "float") then
        kcnt=kcnt+1
        iflag=1
        nsim=nsim+1
        do 5 kk=1,11
5       spar(kk,nsim)=alog10(bexp(k))
c   use a random number generator to dither exponent
        icount=0
4       continue
c        write(89,*)'icount=   X',icount
c        if (icount .ne. 10)  then

          ran=rand()

          idum=int(abs(1000000*ran))
c          ran=5.5*(ran-0.5)
c          ran=2.5*(ran-0.5)

c  if dither is too small, try again
c        if (abs(ran) .lt. 0.2) write(89,*)'ran=',ran,' icount=',icount
c          if (abs(ran) .lt. 0.2)  go to 4
          bexp_orig=bexp(k)
c          write(89,*)'here ',exp_type(k)
          if (exp_type(k) .eq. 'm' ) then
c            write(89,*)'here ',exp_type(k)
            bexp(k)=(alog10(20.0/365.25)-alog10(tmin/(5*365.25)))*ran
     &      + alog10(tmin/(5.*365.25))
          end if
          if (exp_type(k) .eq. 'e' ) then
            bexp(k)=
     &      (alog10(tmax*5./365.25)-alog10(tmin/(5*365.25)))*ran
     &      + alog10(tmin/(5*365.25))
          end if
c          bexp(k)=ran+alog10(bexp(k)) 
c          write(89,*)'ran log(bexp)',ran, bexp(k)
          bexp(k)= 10.0**bexp(k)
c          write(89,*)'  bexp ', bexp(k)
          if ((bexp(k) .lt. tmin/(365.0*10.)) .or. 
     &       (bexp(k) .gt. 365.*10.*tmax))  then
             bexp(k)=bexp_orig
c           write(89,*)'icount = ', icount,'  bexp=',bexp(k)
             icount=icount+1
             go to 4
          end if
c        else
c   got stuck in a loop; just assign random number for bexp()
c         print*, ' stuck in loop; bailed'
c         ran=rand()
c         ran=ran-0.5
         
c         bexp(k)=10**ran
c         write(89,*)' bexp(k)=',bexp(k),' ran=',ran
         

c        end if
c        write(89,*)' testing', k,bexp(k),alog10(bexp(k)),spar(1,nsim)
        if (abs(alog10(bexp(k)) -spar(1,nsim)) .lt. 0.05 ) then
           bexp(k)=bexp(k)*0.12*float(kcnt+1)
           write(89,*)' bexp(k) to close to first one; modify'
           write(89,*) ' bexp(1)= ',bexp(1),' bexp(k)=',k,bexp(k)
        end if

        call modify_A(A,dtime,texp)


        call calcres(0,ModType,sumRes,A,d,res,chi2)

        spar(kcnt+1,nsim)=alog10(bexp(k))
        sres(kcnt+1)=chi2

        bexp(k)=bexp_orig
        nmod=nmod_orig    
        end if

1     continue
      rewind(23)
      write(23,*)idum

c
c  use downhill simplex to get best fit time-constants
c
      ftol=0.000001
      ftol=2.0e-05
c  check to see if we need to iterate to find best fit time constants
      if (iflag .eq. 0) then
        write(13,130)(x(i),e(i),i=1,nmod)
        return
      end if
      do 3 k=1,nsim+1

3       write(89,*) ' input to Amoeba, k=',k,sres(k),
     &       (spar(k,j),10.0**spar(k,j),j=1,nsim)

      call NedlerMeadEXP(spar,sres,11,10,nsim,FTOL,ITER,
     &  ModType,A,d,res,max_data,max_mod)
 
      do 41 k=1,nsim+1

41       write(89,*) ' output from Amoeba, k=',k,sres(k),
     &       (j,spar(k,j),10.0**spar(k,j),j=1,nsim)
         write(89,*)" "


c      if (ITER .ge. 75 ) then
c   Try again !!
c      end if

      fmin=1.e+30
      kcnt=0
      ibest=0
      do 10 i=1,nsim+1

         if (sres(i) .lt. fmin) then
         ibest=i
         fmin=sres(i)
         end if
10    continue

      write(89,*),' ibest=',ibest
      izz=1
      do 11 j=1,nexp
        if (exp_choice(j) .eq. "float" ) then
           bexp(j)=10.0**spar(ibest,izz)
           izz=izz+1
        end if
11    continue
      write(89,*)" Best fit ",fmin,(j, bexp(j),j=1,nexp)
      write(89,*) " "

      nmod=nmod_orig

      call modify_A(A,dtime,texp)
      call calcres(0,ModType,sumRes,A,d,res,chi2)

c   put the estimate of exponential into list of x(i)
         ix=0
         do 18 k=1,nexp
           if (exp_choice(k) .eq. "float") then
           ix=ix+1
           x(nmod+ix)=bexp(k)
           end if
18       continue
c
c  estimate the standard error on all model parameters by taking first differences
c  and doing the appropriate inversions
c  

      if (nopt .eq. 1 ) then     
      do 185 j=1,nmod
185   Anor(j)=1.0


c  next, add columns associated with the time-constants
        kexp=0
        kkexp=0
        do 19 kexp=1,nexp
          if (exp_choice(kexp) .eq. "float") then
          kkexp=kkexp+1

          xxd=1.1
          if (bexp(kexp) .gt. 1.0) xxd=0.9

          Amin=1.0e+20
          Amax=-1.0e+20
          nmod=nmod+1

c            print*,' x(nmod-1)=',x(nmod-1),' nmod-1=',nmod-1, 
c     &  'nmod-nexp+kexp-1=',nmod-nexp+kexp-1,'x()=',x(nmod-nexp+kexp-1)
            do 20 j=1,ic
              A(j,nmod)=0.0
              if (exp_type(kexp) .eq. "e") then
              if (dtime(j) .ge. texp(kexp)-t_start) then

c   use the derrivative here  ( need to be careful about units between days and years for time
c              tuuu=(dtime(j)-(texp(kexp)-t_start))
c              A(j,nmod)=(-1.0*tuuu/( (365.25*bexp(kexp))**2 ))
c     &              * exp(-1.0*tuuu/(365.25*bexp(kexp) ) )
c             A(j,nmod)=A(j,nmod)*x(nmod_orig+kexp) * 365.25
              tuuu=(dtime(j)-(texp(kexp)-t_start))/365.25
              A(j,nmod)=(-1.0*tuuu/( (bexp(kexp))**2 ))
     &              * exp(-1.0*tuuu/(bexp(kexp) ) )
             A(j,nmod)=A(j,nmod)*x(nmod_orig+kkexp)

             end if

             end if
             if (exp_type(kexp) .eq. "o") then
              if (dtime(j) .gt. texp(kexp)-t_start) then


              e1=alog10(abs(bexp(kexp)) +
     &            (sngl(dtime(j)-(texp(kexp)-t_start))/365.25))
              e2=alog10(abs(bexp(kexp))*xxd +
     &            (sngl(dtime(j)-(texp(kexp)-t_start))/365.25))
             A(j,nmod)=1.*(e1-e2)/(1.0*(bexp(kexp)-xxd*bexp(kexp)))
     &                *x(nmod_orig+kkexp)
c             print*,j,bexp(kexp),e1,e2,A(j,nmod)
c             print*,"  ",e1-e2,xxd,x(nmod_orig+kexp)
c             print*,"    ", (bexp(kexp)-xxd*bexp(kexp))
             end if
             end if
             if (exp_type(kexp) .eq. "m" ) then
              if (dtime(j) .gt. texp(kexp)-t_start) then
c              e1=alog10( 1.0 + 
c     &            (sngl(dtime(j)-(texp(kexp)-t_start))/365.25)
c     &             /abs(bexp(kexp)) )
c              e2=alog10( 1.0 + 
c     &            (sngl(dtime(j)-(texp(kexp)-t_start))/365.25)
c     &             /abs(bexp(kexp))*xxd )
c             A(j,nmod)=1.*(e1-e2)/(1.0*(bexp(kexp)-xxd*bexp(kexp)))
c     &                *x(nmod_orig+kkexp)
c             print*,j,nmod,A(j,nmod),x(nmod_orig+kkexp),bexp(kexp)
             tuuu=(dtime(j)-(texp(kexp)-t_start))/365.25
             A(j,nmod)=x(nmod_orig+kkexp)*(-1.0)*tuuu/
     &          (((bexp(kexp))**2)*(1.0+tuuu/bexp(kexp))*alog(10.0))
c             print*,tuuu,A(j,nmod)

             end if
             end if
             if (A(j,nmod) .lt. Amin) Amin=A(j,nmod)
             if (A(j,nmod) .gt. Amax) Amax=A(j,nmod)
c          write(87,*)j,t(j),(A(j,k),k=1,nmod),bexp(kexp),e1,e2
c          write(86,*)(A(j,k),k=1,nmod)

20         continue
           range=Amax-Amin
           Anor(nmod)=200.0/range
c           Anor(nmod)=1.0/abs(x(nmod_orig+kexp))
c           Anor(nmod)=1.0/(range*abs(x(nmod_orig+kexp)))
           Anor(nmod)=1.0/(range)
c           print*,nmod,kexp,' nmod and kexp'
           print*,"max, min, range", Amax,Amin,range,' Norm',
c     &    Anor(nmod),' x(nmod-2*nexp+kkexp)=',x(nmod-2*nexp+kkexp),
c     &    " (nmod-2*nexp+kkexp)= ",(nmod-2*nexp+kkexp),nmod,2*nexp,
c     &     kkexp,
     &    Anor(nmod),' x(nmod_orig+kkexp)=',x(nmod_orig+kkexp),
     &    " (nmod_orig+kkexp)=",nmod_orig+kkexp
c   renormalize the nmod column
           do 21 j=1,ic
21         A(j,nmod)=A(j,nmod)*Anor(nmod)
c     if bexp is 'out of range' zero out the A matrix

           if ((bexp(kexp) .le. tmin/(365.0*10.)) .or.
     &         (bexp(kexp) .gt. 5.*tmax/365.))  then
               do 213 j=1,ic
213            A(j,nmod)=0.0
           end if
          end if
19      continue 
      mdmax=max_data
      maxnmod=max_parm
c
c  Compute model covariance depending upon ModType
c
      if ( ModType .eq. 'f' ) then
c        print*,' I am here', (Finv(i),i=1,10)
c        do 890 i=1,10
c        print*,i,(B(i,j),j=1,10)
c890     continue
c     whiten model
        do 71 j=1,nmod
          do 72 i=1,ic
            din(i)=(A(i,j))
72        continue
          call convolv1(ic,Finv,din,dout,max_time,mdmax)
          do 73 i=1,ic
            Aw(i,j)=(dout(i))
            Awt(j,i)=Aw(i,j)
73        continue

71      continue
c
c    Aw^t * Aw
c

         maxnmod=max_mod
         call dgemm('T','N',nmod,nmod,ic,1.0d+0,Aw,mdmax,Aw,mdmax,
     & 0.0d+0,a2,maxnmod)

         if (nmiss .gt. 0) then
c    Aw^t B Aw
c    a3 = Aw^t B
           call dgemm('T','N',nmod,ic,ic,1.0d+0,Aw,mdmax,B,mdmax,
     &     0.0d+0,a3,maxnmod)
c a2= (Aw^t B) Aw= a2 - a3 Aw
           call dgemm('N','N',nmod,nmod,ic,-1.0d+0,a3,maxnmod,Aw,mdmax,
     &   1.0d+0,a2,maxnmod)


         end if
c
c invert Aw^t * Aw
c

 
c  end covariance for ModType f
      else
c
c  Covariance for Modtype c and n (Cholesky decomp).
c
c  compute At*covinv
c
      call dgemm('T','N',nmod,ic,ic,1.0d+0,A,max_data,covinv,max_data,
     & 0.0d+0,a1,max_mod)
c
c compute At*covinv*A   A1*A
c
      call dgemm('N','N',nmod,nmod,ic,1.0d+0,a1,max_mod,A,max_data,
     & 0.0d+0,a2,max_mod)
c
c  End switch for ModType
c
      end if
      do 120 i=1,nmod
      do 119 j=1,nmod
119    ad2(i,j)=a2(i,j)
120    continue
c
c invert At*covinv*A
c
      call dsyev('V','L',nmod,ad2,max_mod,evd,wzd,5000,ier)

      print*,"In time constant routine,  ier =",ier
      print*,' eigenvalues',(evd(k),k=1,nmod)
      ixx=0
      do 2221 k=nmod,1,-1
        if (abs(evd(k)/evd(nmod)) .gt. 1.0e-09) ixx=ixx+1
2221  continue
      do 220 i=1,nmod

      do 221 j=1,nmod
        at=0.0
        do 222 k=nmod,nmod-ixx+1,-1
c         if (abs(evd(k)/e(nmod)) .gt. 1.0e-07) then
            at=at+ad2(i,k)*ad2(j,k)/evd(k)
c            ixx=ixx+1
c         end if
c        if (e(k) .gt. 0) 
c    &     at=at+a2(i,k)*a2(j,k)/e(k)
222     continue
221   ai(i,j)=at
c      print*," ai (i,j) i=",i,(ai(i,j),j=1,nmod)
220   continue
      print*,' Using',ixx,' out of', nmod,' eigenvalues'


      do 330 i=1,nmod
c      print*,'i=',i,' ai=',ai(i,i),' sqrt',sqrt(ai(i,i)),
c     &   ' Anor=',Anor(i)
      write(6,2226) i,ai(i,i), sqrt(ai(i,i)), Anor(i)
2226  format (' i=',i5,'  ai= ',e12.5,'  sqrt=',e12.5,'  Anor=',e10.3)
c     &,sqrt(ai(i,i))*Anor(i)
      if (ai(i,i) .gt. 0) then
         e(i)=sqrt(ai(i,i))*Anor(i)
         else
         e(i)=1.0e+06
c         e(i)=sqrt(-ai(i,i))*Anor(i)
      end if
c      print*,sqrt(ai(i,i))*Anor(i),e(i)
330   continue
c
c  Output the model Covariance
c
      rewind (78)
      write(78,319)
319   format("  Model Covariance",/," Listed approximately in same",
     & "   order as output with exception of floating ",
     & /,"   omori/exponental amplitude and time constant",/)
      do 320 i=1,nmod
        write(78,340)i,sqrt(ai(i,i))*Anor(i),
     &  (ai(i,j)*((Anor(i)*Anor(j))**2),j=1,nmod)
    
320   continue
340   format(i3,' std. err= ',e13.4,' cov= ',82e13.4)
c  Output correlation coefficients
      write(78,318)
318   format(//,"  Model correlation coef.",/," Listed approximately ",
     & "in same",
     & "   order as output with exception of floating ",
     & /,"   omori/exponental amplitude and time constant",/)     
      do 321 i=1,nmod
        write(78,341)i,sqrt(ai(i,i))*Anor(i),
     &  (ai(i,j)/(sqrt(ai(i,i))*sqrt(ai(j,j))),j=1,nmod)
    
321   continue
341   format(i3,' std. err= ',e13.4,' Cross correl ',82f8.4)
      end if
c
c  END of error calculation
c
c      if ( iswitch .eq. 1) then
      if ( nopt .eq. 1) then
        write(13,130)(x(i),e(i),i=1,nmod)


c        write(88,*)" ERROR calc"
c130     format(90(3x,f15.3,1x,f15.3))
130     format(90(3x,e17.10,1x,e16.10))
      else
        write(13,130)(x(i),e(i),i=1,nmod),(x(nmod+i),0.0,i=1,ix)
c        write(88,*)" no error"
      end if

      nmod=nmod_orig
      return
      end
      subroutine calcres(iswitch,ModType,sum,ax,d,dr,chi2)

c max_data  dimension of data
c max_parm  dimension of num parameters
c ic number of data
c nmod number of parameters
c ax  design matrix
c d  data
c dr  residuals
c covinv  inverse covariance matrix
c covar  covariance matrix
c iswitch  outputs the model and its error if set to 1
c  x  the model
c  e  the error
c  sum Sum of residuals normalized by covariance  dr^t*covinv*dr
c  ModType;  f, c, n
c

      double precision ax(max_data,max_mod),dtime(12687),
     &  covinv(12687,12687),
     & covar(12687,12687),texp(82),t_start,t_stop,a1(max_mod,max_data),
     & a2(max_mod,max_mod),wzd(5000),evd(max_mod),t_year(12687),
     &  at,ai(82,82),ainv(82,12687),chi2,resA,Finv(32768)
      double precision Aw(max_data,max_mod),B(12687,12687),
     & a3(max_mod,max_data),Awt(max_mod,max_data),sum
      dimension d(max_data),dr(max_data),x(82),e(82),bexp(82)
      double precision dw(max_data),din(max_data),dout(max_data)
      integer irowmiss(12687)
      character*1, ModType
c   double precision
      common /ModFit1/dtime,t_year,covinv,covar,texp
      common /ModFit1a/t_start,t_stop
c  single precision
      common /ModFit2/x,e,bexp
      common /ModFit2a/expmax
c integer
      common /ModFit3/max_data,max_mod,ic,nmod,n_exp,nmiss,irowmiss,
     & max_time,ipl_flag_1,ipl_flag_2,ibp_flag,nmiss_max

c  character
c      common /ModFit4/ exp_choice,exp_type
c  double precision
      common /Modfit5/ Finv, B
      mdmax=max_data
      maxnmod=max_mod
c      print*,' In calcres',iswitch,t_start,t_stop
c      print*,max_data,max_mod
c      print*,(i,covinv(i,i),i=1,10)
c      print*,(i,covinv(i,i),i=ic-10,ic+2)
c      print*,(ax(1,j),j=1,12)
c      print*,(ax(ic,j),j=1,12)
c-------------------------------------------------------------
c
c  find new model that fits data using updated covariance matrix
c


c      if (( ModType .eq. 'f' ) .and. (iswitch .ne.1) ) then
      if (( ModType .eq. 'f' ) ) then
c
c      execute this by "whitening" both the model and data using Finv
c
c       C^-1 = (F^-1)^t * (F^-1)
c       x=A^-1 d = (A^t C^-1 A)^-1 A^t C^-1 d
c                = [ A^t (F^-1)^t * (F^-1) A ]^-1 A^t (F^-1)^t (F^-1 d)
c       let F^-1 d = dw
c           and (F^-1) A = Aw

c
c      x= [Aw^t Aw - Aw^t B Aw ]^-1 [Aw^t - Aw^B] dw
c   where B= E (E ^t E)^-1 E^t
c   and E = F^-1 M  where M is a matrix composed of indices of missing data
c   Note the B is computed in funmin and passed onto 
c
c   chi^2 =d^t * C^-1 d
c         =dw^t dw -dw^t *B dw

c
c
c     whiten model
        do 71 j=1,nmod
          do 72 i=1,ic
            din(i)=(Ax(i,j))
72        continue
          call convolv1(ic,Finv,din,dout,max_time,mdmax)
          do 73 i=1,ic
            Aw(i,j)=(dout(i))
            Awt(j,i)=Aw(i,j)
73        continue

71      continue
         
c     whiten data
        do 711 i=1,ic
711     din(i)=dble(d(i))
        call convolv1(ic,Finv,din,dw,max_time,mdmax)


c
c    Aw^t * Aw
c


         call dgemm('T','N',nmod,nmod,ic,1.0d+0,Aw,mdmax,Aw,mdmax,
     & 0.0d+0,a2,maxnmod)


      if (nmiss .gt. 0) then
c    Aw^t B Aw
c    a3 = Aw^t B
         call dgemm('T','N',nmod,ic,ic,1.0d+0,Aw,mdmax,B,mdmax,
     &   0.0d+0,a3,maxnmod)


c a2= (Aw^t B) Aw= a2 - a3 Aw
         call dgemm('N','N',nmod,nmod,ic,-1.0d+0,a3,maxnmod,Aw,mdmax,
     &   1.0d+0,a2,maxnmod)

c    Awt = Aw^t - Aw^t B = Aw^t - a3
       do 666 i=1,nmod
       do 667 j=1,ic
         Awt(i,j)=Awt(i,j)-a3(i,j)
667    continue
666    continue

      end if
c
c invert Aw^t * Aw
c

        call dsyev('V','L',nmod,a2,maxnmod,evd,wzd,5000,ier)


        if (iswitch .eq. 1) print*,' Eigenvalues ',(i,evd(i),i=1,nmod)
c examine the eigenvalues
        ie=1
        iflag=0
        do 212 i=nmod-1,1,-1
          if (iflag .eq. 1) go to 212
          rat=evd(i)/evd(nmod)
          rat=abs(rat)
          if (rat .lt. 1.0e-09) go to 213
          ie=ie+1
          go to 212
213       iflag =1
212     continue
        if (iswitch .eq. 1)
     &    print*,' Using ',ie,' out of', nmod,' eigenvalues'
        if (ie .ne. nmod) 
     &    print*,' Using ',ie,' out of', nmod,' eigenvalues'
c  construct inverse
     
        if (ier .ne. 0) then

          print*,' eigenvalues ',(evd(i),i=1,nmod)
          print*,' ier= ',ier
          print*,'  if ier=0; success'
          print*,'        <0the i-th argument had an illegal value'
          print*,'        >0the algorithm  failed  to  converge'
          print*,'      do a man ssyev'
        end if
        do 223 i=1,nmod
        do 224 j=1,nmod
          at=0.0
          do 225 k=nmod,nmod-ie+1,-1
225       at=at+a2(i,k)*a2(j,k)/(1.0*evd(k))
224     ai(i,j)=at
223     continue


c
c compute the inverse matrx 
c

        call dgemm('N','N',nmod,ic,nmod,1.0d+0,ai,maxnmod,Awt,maxnmod,
     &   0.0d+0,ainv,maxnmod)

        do 430 i=1,nmod

        x(i)=0.
        e(i)=0.
           do 431 j=1,ic
431        x(i)=x(i)+ainv(i,j)*dw(j)

430     continue

c
c compute the residuals
c
        do 335 i=1,ic
        dr(i)=0.
          do 336 j=1,nmod
336       dr(i)=dr(i)+ax(i,j)*x(j)

        din(i)=dble(dr(i))
        dr(i)=d(i)-dr(i)
335     continue


c     whiten residuals

        call convolv1(ic,Finv,din,dout,max_time,mdmax)
        chi2=0.0
        do 76 i=1,ic
          chi2=chi2+(dout(i)-dw(i))**2

76      continue
c
c      chi2adj= (dout-dw)^t * B * (dout-dw)
c  
        do 7687 i=1,ic
          din(i)=0.0
          do 7686 j=1,ic
          din(i)=din(i)+B(i,j)*(dout(j)-dw(j))
7686     continue
7687    continue
        chi2adj=0.0
        do 7688 i=1,ic

          chi2adj=chi2adj+ (dout(i)-dw(i))*din(i)

7688    continue



        chi2=chi2-chi2adj         
c     
      else

c   execute the following using standard inversion methods
c
c
c  compute At*covinv
c

      call dgemm('T','N',nmod,ic,ic,1.0d+0,ax,mdmax,covinv,mdmax,
     & 0.0d+0,a1,maxnmod)

c
c compute At*covinv*A   A1*A
c
      call  dgemm('N','N',nmod,nmod,ic,1.0d+0,a1,maxnmod,ax,mdmax,
     & 0.0d+0,a2,maxnmod)


c
c invert At*covinv*A
c

      call dsyev('V','L',nmod,a2,maxnmod,evd,wzd,5000,ier)

      if (iswitch .eq. 1) print*,' Eigenvalues ',(i,evd(i),i=1,nmod)
c examine the eigenvalues
      ie=1
      iflag=0
      do 210 i=nmod-1,1,-1
        if (iflag .eq. 1) go to 210
        rat=evd(i)/evd(nmod)
        rat=abs(rat)
        if (rat .lt. 1.0e-09) go to 211
        ie=ie+1
        go to 210
211     iflag =1
210   continue
       if (iswitch .eq. 1)
     &   print*,' Using ',ie,' out of', nmod,' eigenvalues'
c  construct inverse
     
      if (ier .ne. 0) then

        print*,' eigenvalues ',(evd(i),i=1,nmod)
        print*,' ier= ',ier
        print*,'  if ier=0; success'
        print*,'        <0the i-th argument had an illegal value'
        print*,'        >0the algorithm  failed  to  converge'
        print*,'      do a man ssyev'
      end if
      do 220 i=1,nmod
      do 221 j=1,nmod
        at=0.0
        do 222 k=nmod,nmod-ie+1,-1
222     at=at+a2(i,k)*a2(j,k)/evd(k)
221   ai(i,j)=at
220   continue

c
c compute the inverse matrx 
c
      call dgemm('N','N',nmod,ic,nmod,1.0d+0,ai,maxnmod,a1,maxnmod,
     & 0.0d+0,ainv,maxnmod)

      do 230 i=1,nmod

      x(i)=0.
      e(i)=0.
         do 231 j=1,ic
231      x(i)=x(i)+ainv(i,j)*d(j)
230   continue
c  output weights that are used to estimate the various parameters --
c   experimental
c      open(66,file='wght.dat')
c      rewind (66)
c      do i=1,ic
c        write(66,*)i,(ainv(k,i),k=1,nmod)
c      end do
c      close (66)
c      print*,' weight matrix in wght.dat'

c
c compute the residuals
c
      do 235 i=1,ic
      dr(i)=0.
        do 236 j=1,nmod
236     dr(i)=dr(i)+ax(i,j)*x(j)
      dr(i)=d(i)-dr(i)
235   continue

c
c  calculate RMS residuals
c

      chi2=0.
      do 3 i=1,ic
      resA=0.
        do 4 j=1,ic
        resA=resA+dble(dr(i))*covinv(i,j)*dble(dr(j))
4       continue
c      print*,i,res(i)
      chi2=chi2+resA
3     continue


      end if
c ------------End of Inversion
c
c
c compute the error bars
c
     

      if (iswitch .eq. 1) then



      do 330 i=1,nmod
c
c   model cov = [A^t * C^-1 A ]^-1  which is ai from above
c

      e(i)=sqrt(ai(i,i))

330   continue
c
c  Output the model Covariance
c
      rewind (78)
      write(78,319)
319   format("  Model Covariance",/," Listed approximately in same",
     & "   order as output with exception of floating ",
     & /,"   omori/exponental amplitude and time constant",/)
      do 320 i=1,nmod
        write(78,340)i,sqrt(ai(i,i)),
     &  (ai(i,j),j=1,nmod)
    
320   continue
340   format(i3,' std. err= ',e13.4,' cov= ',82e13.4)
c  Output correlation coefficients
      write(78,318)
318   format(//,"  Model correlation coef.",/," Listed approximately ",
     & "in same",
     & "   order as output with exception of floating ",
     & /,"   omori/exponental amplitude and time constant",/)     
      do 321 i=1,nmod
        write(78,341)i,sqrt(ai(i,i)),
     &  (ai(i,j)/(sqrt(ai(i,i))*sqrt(ai(j,j))),j=1,nmod)
    
321   continue
341   format(i3,' std. err= ',e13.4,' Cross correl ',82f8.4)

      end if

      return
      end
      subroutine output(x,e,max_mod,net,nbt,nrate,n_rat_chg,
     & nper,noff,rat_chg1,rat_chg2,off,per,sea,rate_norm,
     &  n_file_press,aux_norm,n_exp,texp,bexp,exp_choice,exp_type,
     & rat_chng_norm)
c  output results of model
      dimension x(max_mod),e(max_mod)
      dimension bexp(10)
      character*7 exp_choice(10),exp_type(10)
      character*3 net
      double precision rat_chg1(82),rat_chg2(82),off(82),per(82),
     &   texp(10), rat_chng_norm(82),dec_time1,dec_time2,
     &  aux_norm(max_mod)
c
c  five different output format; otr, otd, otx,mjd, and gmt
c
      nc=0
      sea=0.

      do 64 i=1,nbt
        write(6,201)i,x(i),e(i)
201     format(' Nomimal value for baseline ',i3,f10.2,' +/- ',f10.2)
64    continue
      nc=nc+nbt
      if (nrate .eq. 1) 
     &      write(6,202)x(nc+1)/rate_norm,e(nc+1)/rate_norm
202   format(' Rate in units per year ',f16.4,' +/- ',f16.4)
      nc=nc+nrate
      if (n_rat_chg .gt. 0) then
        do 65 k=1,n_rat_chg
        itime1=int(rat_chg1(k))
        dec_time1=rat_chg1(k)-float(itime1)
        itime2=int(rat_chg2(k))
        dec_time2=rat_chg2(k)-float(itime2)

        if (net .eq. 'otr') then
          call inv_jul_time(itime1,nyr,jul)
          call inv_jul_time(itime2,nyr2,jul2)
          write(6,203)k,nyr,jul+dec_time1,nyr2,jul2+dec_time2,
     &       x(nc+k)/rat_chng_norm(k),e(nc+k)/rat_chng_norm(k)
        end if
        if (net .eq. 'otx' ) then
          call invcal(itime1,nyr,mn,idate)
          call invcal(itime2,nyr2,mn2,idate2)
          write(6,204)k,nyr,mn,idate+dec_time1,
     &      nyr2,mn2,idate2+dec_time2,
     &       x(nc+k)/rat_chng_norm(k),e(nc+k)/rat_chng_norm(k)
        end if
        if (net .eq. 'otd') then
          call invcal(itime1,nyr,mn,idate)
          call invcal(itime2,nyr2,mn2,idate2)
          write(6,2041)k,nyr,mn,idate,dec_time1,
     &      nyr2,mn2,idate2,dec_time2,
     &       x(nc+k)/rat_chng_norm(k),e(nc+k)/rat_chng_norm(k)
        end if
        if (net .eq. 'gmt') then

          call invcal(itime1,nyr,mn,idate)
          call invcal(itime2,nyr2,mn2,idate2)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          ihr2=int(24.0*dec_time2)
          imn2=int(1440.d+0*(dec_time2-dble(ihr2)/24.0d+0))
          sec2=dec_time2-dble(ihr2)/24.0d+0-dble(imn2)/(1440.0d+0)
          sec2=sec2*3600.0*24.0
          isec2=int(sec2)
          print*,sec2,isec2
          sec2=sec2-float(isec2)
          write(6,2042)k,nyr,mn,idate,ihr1,imn1,isec1,sec1,
     &      nyr2,mn2,idate2,ihr2,imn2,isec2,sec2,
     &       x(nc+k)/rat_chng_norm(k),e(nc+k)/rat_chng_norm(k)
        end if    
        if (net .eq. 'mjd' ) then
          write(6,2043)k,dble(itime1)+dec_time1+36933.d+0,
     &     dble(itime2)+dec_time2+36933.d+0, 
     &       x(nc+k)/rat_chng_norm(k),e(nc+k)/rat_chng_norm(k)

        end if
65      continue
      end if
      nc=nc+n_rat_chg
203   format(' Rate change number ',i3,' between ',i5,f10.3,
     & ' and ',i5,f10.3,
     & ' is ',f16.4,' +/- ',f16.4)
204   format(' Rate change number ',i3,' between ',2i5,f6.3,
     & ' and ',2i5,f6.3,
     & ' is ',f12.4,' +/- ',f10.4)
2041   format(' Rate change number ',i3,' between ',i4,i2.2,i2.2,f4.3,
     & ' and ',i4,i2.2,i2.2,f4.3,
     & ' is ',f12.4,' +/- ',f10.4)
2042   format(' Rate change number ',i3,' between ',i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1,
     & ' and ',i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1,
     & ' is ',f12.4,' +/- ',f10.4)
2043   format(' Rate change number ',i3,' between ',f16.9,' and ',
     & f16.9,' is ',f12.4,' +/- ',f10.4)
      if (nper .gt. 0) then
         do 66 k=1,nper
         fmag=x(nc+2*k-1)**2 + x(nc+2*k-0)**2
         fmag=sqrt(fmag)
         e_mag=(x(nc+2*k-1)**2)*(e(nc+2*k-1)**2)
     &    +(x(nc+2*k-0)**2)*(e(nc+2*k-0)**2)
         e_mag=sqrt(e_mag)/fmag
         write(6,205) per(k),x(nc+2*k-1),e(nc+2*k-1), 
     &   x(nc+2*k-0),e(nc+2*k-0), fmag,e_mag
         if (k .eq. 1) sea=fmag
66       continue
      end if
      nc=nc+2*nper
205   format(' Period of ',f8.3,' days,  cos amp= ',f12.2,' +/-',f10.2,
     & '  sin amp= ',f12.2,' +/-',f10.2,'  magnitude= ',
     &  f12.2,' +/-',f10.2)
      if (noff .gt. 0) then
        do 67 k=1,noff
        itime=int(off(k))
        dec_time1=off(k)-float(itime)
        if (net .eq. 'otr') then
          call inv_jul_time(itime,nyr,jul)
          write(6,206)k,nyr,jul+dec_time1,x(nc+k),e(nc+k)
        end if
        if (net .eq. 'otx') then
          call invcal(itime,nyr,mn,idate)
          write(6,207)k,nyr,mn,idate+dec_time1,x(nc+k),e(nc+k)
        end if
        if (net .eq. 'otd') then
          call invcal(itime,nyr,mn,idate)
          write(6,2071)k,nyr,mn,idate,dec_time1,x(nc+k),e(nc+k)
         end if
        if (net .eq. 'gmt') then
          call invcal(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,2072)k,nyr,mn,idate,ihr1,imn1,isec1,sec1,
     &       x(nc+k),e(nc+k)
        end if
        if (net .eq. 'mjd') then
          write(6,2073)k,dble(itime)+dec_time1+36933.d+0,
     &       x(nc+k),e(nc+k)

        end if
67      continue
      end if
206   format(' Offset number ',i5,' at ',i5,f10.3,
     & ' is ',f10.2,' +/- ',f8.2)
207   format(' Offset number ',i3,' at ',i5,i3,f8.3,
     & ' is ',f10.2,' +/- ' ,f8.2)
2071  format(' Offset number ',i3,' at ',i4,i2.2,i2.2,f4.3,
     & ' is ',f10.2,' +/- ' ,f8.2)
2072   format(' Offset number ',i3,' at ',i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1,
     & ' is ',f12.4,' +/- ',f10.4)
2073   format(' Offset number ',i3,' at ',f16.9,
     & f12.4,' +/- ',f10.4)
      nc=nc+noff
      ikk=n_exp
c      nc=nc+n_file_press
      nfix=0      
      if (n_exp .ne. 0) then
c      do 680 k=1,n_exp
c680   if (exp_choice(k) .eq. "float") nfloat=nfloat+1
c      print*,"nfloat=",nfloat
      do 68 k=1,n_exp

         if (exp_choice(k) .ne. "float") then
         
         nfix=nfix+1
         itime=int(texp(k))
         dec_time=texp(k)-float(itime)
        if (net .eq. 'otr') then
          call inv_jul_time(itime,nyr,jul)
          write(6,208)k,nyr,jul+dec_time,x(nc+k),e(nc+k),bexp(k)
        end if
        if (net .eq. 'otx') then
          call invcal(itime,nyr,mn,idate)
          write(6,209)k,nyr,mn,idate+dec_time,x(nc+k),e(nc+k),bexp(k)
        end if
        if (net .eq. 'otd') then
          call invcal(itime,nyr,mn,idate)
          write(6,210)k,nyr,mn,idate,dec_time,x(nc+k),e(nc+k),bexp(k)
        end if
        if (net .eq. 'gmt') then
          call invcal(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,211)k,nyr,mn,idate,ihr1,imn1,isec1,sec1,
     &      x(nc+k),e(nc+k),bexp(k)
        end if
        if (net .eq. 'mjd') then
           write(6,212)k,dble(itime)+dec_time1+36933.d+0,
     &       x(nc+k),e(nc+k),bexp(k)
        end if
        else
         ikk=ikk+1
         nc=nc+n_file_press
         itime=int(texp(k))
         dec_time=texp(k)-float(itime)
        if (net .eq. 'otr') then
          call inv_jul_time(itime,nyr,jul)
          write(6,2081)k,nyr,jul+dec_time,x(nc+k),e(nc+k),
     &      bexp(k),e(nc+0+ikk)
        end if
        if (net .eq. 'otx') then
          call invcal(itime,nyr,mn,idate)
          write(6,2091)k,nyr,mn,idate+dec_time,x(nc+k),e(nc+k),
     &      bexp(k),e(nc+0+ikk)
        end if
        if (net .eq. 'otd') then
          call invcal(itime,nyr,mn,idate)
          write(6,2101)k,nyr,mn,idate,dec_time,x(nc+k),e(nc+k),
     &      bexp(k),e(nc+0+ikk)
        end if
        if (net .eq. 'gmt') then
          call invcal(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,2111)k,nyr,mn,idate,ihr1,imn1,isec1,sec1,
     &      x(nc+k),e(nc+k),bexp(k),e(nc+0+ikk)
        end if
        if (net .eq. 'mjd') then
           write(6,2121)k,dble(itime)+dec_time1+36933.d+0,
     &       x(nc+k),e(nc+k),bexp(k),e(nc+0+ikk)
        end if
         nc=nc-n_file_press

        end if
68      continue 
208   format(' Exponential number ',i5,' at ',i5,f10.3,
     & ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ', e12.5,' yrs')
c     & ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ',f10.5,' yrs')

209   format(' Exponential number ',i3,' at ',i5,i3,f8.3,
     & ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ', e12.5,' yrs')
210   format(' Exponential number ',i3,' at ',i4,i2.2,i2.2,f4.3,
     & ' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs')  
211   format(' Exponential number ',i3,' at ',i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1,
     & ' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs')
212   format(' Exponential number ',i3,' at ',f16.9,
     & f12.4,' +/- ',f10.4,'  Time constant is: ', e12.5,' yrs')
2081  format(' Exponential number ',i5,' at ',i5,f10.3,
     & ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ',e10.5,
     &  ' +/- ',e11.5,' yrs')
2091  format(' Exponential number ',i3,' at ',i5,i3,f8.3,
     & ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ',e10.5,
     &  ' +/- ',e11.5,' yrs')    
2101  format(' Exponential number ',i3,' at ',i4,i2.2,i2.2,f4.3,
     & ' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs',
     &  ' +/- ',e11.5,' yrs') 
2111  format(' Exponential number ',i3,' at ',i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1,
     & ' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs',
     &  ' +/- ',e11.5,' yrs')     
2121  format(' Exponential number ',i3,' at ',f16.9,
     & f12.4,' +/- ',f10.4,'  Time constant is: ', e12.5,' yrs',
     &  ' +/- ',e11.5,' yrs') 
      end if
c      print*,'n_file_press=',n_file_press,' nc=',nc
c      nc=nc-n_file_press
c  don't need to add to nc; the exponental is 'out of order'; nc was added to in noff.
c      nc=nc
c      nc=nc+n_exp-ikk
c      print*,' nc=',nc,' ikk=',ikk
c      nc=nc-ikk
c      print*,' new nc after adding n_exp', nc
      nc=nc+nfix
      if (n_file_press .ne. 0) then
        do 679 k=1,n_file_press

c        print*,'nc+k',nc+k,'x(nc+k)',x(nc+k),'aux_norm(k)',aux_norm(k)
        write(6,2108)k,x(nc+k)/aux_norm(k),e(nc+k)/aux_norm(k)
c        print*,k,x(nc+k),e(nc+k),aux_norm(k)
679     continue
      end if
2108  format(' response of ', i3, ' input function',f13.5,' +/- ',f12.5)


      return
      end

      subroutine outResid(ic,nmod,A,d,x,dtime,t_start,net,ndim,mdim)

c
c  Calculates and Outputs residuals in different formats using flag net
c
c  A -- design matrix (dimensioned ndim by mdim)
c  d --- data vector
c  x -- the model parameters
c  dtime == time, day number
      dimension d(ndim),x(mdim)
      double precision dtime(ndim),t_start,dec_timed,A(ndim,mdim)
      character*3 net
c      print*,' in outResid'
c      print*,(i,x(i),i=1,nmod)
      print*,'ic=', ic

      do 9011 i=1,ic
        calc_m=0.
        do 9012 j=1,nmod
        calc_m=calc_m+A(i,j)*x(j)
9012    continue
      res_m=d(i)-calc_m
      itime=int(dtime(i)+t_start)
      dec_timed=(dtime(i)+t_start)-int(dtime(i)+t_start)+1.1564e-07
      if (net .eq. 'otr' ) then
         call inv_jul_time(itime,nyr,jul)

          write(65,6501)float(nyr),float(jul)+dec_timed,
     &   res_m,calc_m,d(i)
      end if
      if (net .eq. 'otx' ) then
        call invcal(itime,nyr,mn,idate)
        write(65,6502)nyr,mn,idate+dec_timed,
     &   res_m,calc_m,d(i)
      end if
      if (net .eq. 'otd') then
        call invcal(itime,nyr,mn,idate)
        write(65,6503)nyr,mn,idate,dec_timed,
     &   res_m,calc_m,d(i)
      end if
      if (net .eq. 'gmt') then
        call invcal(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_timed)
          imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
          sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-isec1
          write(65,6504)nyr,mn,idate,ihr1,imn1,isec1,sec1,
     &   res_m,calc_m,d(i)
      end if
      if (net .eq. 'mjd' ) then
         write(65,6505)dble(itime)+dec_timed+36933.d+0,
     &   res_m,calc_m,d(i)
      end if

9011  continue
6501  format(1x,f6.0,1x, f15.9, 3(1x,f15.2))
6502  format(1x,2i5,f15.9, 3(1x,f15.2))
6503  format(1x,i4,i2.2,i2.2,f10.9, 3(1x,f15.2))
6504  format(1x,i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1, 3(1x,f15.2))
6505  format(1x,f18.9, 3(1x,f15.2))
      return
      end
