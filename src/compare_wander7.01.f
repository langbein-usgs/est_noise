c compares the RMS wander or drift of real data with that from from synthetic data
c   which is generated as either
c  uses two power law functions, white noise, and band-pass filtered white noise
c
c  USED IN COMPANION with est_noise7.2x
c
c
c  1) program 1st computes the 'wander' for the real data more or less defined
c      by Agnew, (1992?) as function of averaging interval
c      averaging interval
c  2)  then loops through as series of synthetic data having the
c       same sampling scheme as the real data
c       computes both wander  of the synthetic data
c  3) outputs the wander for the real data, the statistic from synthetic data.
c       The statists are the confidence interval (68, 90, and 95) based upon ranking the
c        results of simulation. In addition, the chi^2 confidence is estimated

c  Modified April 2004 for differences between est_noise5fd and est_noise5 in computing
c    covariances.
c
c Modify Nov 2013 for double precision time
c
c  This is a re-write of compare_wander6  incorporating two different methods to compute
c  the temporal correlations; either by adding several independent noise sources, as
c   in the original est_noise6/compare_wander6, ie quadrature, or adding the filter functions
c    representing different noise and convolving with a single random noise generator,
c    known as 'additive noise' compatable with est_noise7.2x with 'a'.
      dimension d(12687), tim_wan(150),wan_dat(150),H(131072),
     & ran1(131072),ran(131072),dat(131072),wan_sim(150)
      dimension save_wan(501,150),temp(501)
      double precision pl1(131072),pl2(131072),wh(131072),bp(131072),
     &  t(12687),timex,t0
      double precision dt_sam
      character*20 filename,name
      character*132 string
      character*3 net
      character*1 ModType
      common /NoiseMod/exp1,amp1,alpha,exp2,amp2,sig,
     & fl,fh,np,amp_bp,fmiss
      common /Filter/wh,pl1,pl2,bp
      open(25,file='wand_in.jrn')
      open(23,file='seed.dat')
      read(23,*)iseed
      max_wan=150
      max_flt=131072
      max_time=max_flt
      max_data=12687
      max_loop=501
c missing data
      fmiss=99999.
      print*,'Input the data type for  processing'

      print*,' otr=data'
      print*,'       format of year, day of year, res, calc, obs'
      print*,' otd=data'
      print*,'       format of YearMnDa,  res, calc, obs'
      print*,' otx=data'
      print*,'       format of year mo da, res, calc, obs'
      print*,' mjd=data'
      print*,'       Modified Julian day, res, calc, obs'
      print*,'For all of the above, the da is double precision'
      print*,'    which allows decimal days to accomodate seconds'
      print*,' '
      print*,' gmt=data with GMT format for time'
      print*,'        year-mo-daThr:mn:secs.x , res, calc, obs'
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
      print*,' input the file name of data residuals'
      print*,'  format of data yr jul_day data; 1999. 80.345 8.1'
      read(5,101) filename
101   format(a20)
      open(21,file=filename)
      write(25,101)filename
      ic=0
1     continue
        read(21,102,end=2)string
102     format(a132)
        call GetData(net,string,timex,dist,err)
        ic=ic+1
        if (ic .eq. 1) t0=timex
        t(ic)=timex-t0
        d(ic)=dist
        go to 1
2     continue
      close(21)
      print*,' Number of observations read',ic
      if (ic .gt. max_data) then
        print*,' Number of data exceeds array dimension of',max_data
        stop
      end if
      print*,' input name of file with noise parameters'
      read(5,101)filename
      open(22,file=filename)
      read(22,7879)name,fmax,idec,iccc,nmod,
     &  sig,
     &  amp1,    
     &  exp1,     
     &  alpha,    
     &  amp_bp,    
     &  amp2,    
     &  exp2
7879  format(a20,3x,f10.3,3i6,7(3x,f10.3,1x,10x))
      close(22) 
      print*,' Noise parameters read from file'
      print*,sig,amp1,exp1,alpha,amp_bp,amp2,exp2
c  get smallest sampling interval
      dt_sam=999999.
      do 3 i=2,ic
3     if (t(i)-t(i-1) .lt. dt_sam) dt_sam=t(i)-t(i-1)
      print*,' Sampling interval is ',dt_sam,' days'
c
c  calculate the so called "wander" of the data for
c  specified periods
c
c  determine the periods; use constant log(time) to
c   increment the periods
c
      tmin=dt_sam
c      tmin=0.25
      tmax=(t(ic)-t(1))
      dtt=0.05
4     continue
      dtt1=1/dtt
c      nt=int(dtt1*(log10(tmax)-log10(tmin) )) + 1
      nt=int(dtt1*(log10(tmax)-log10(2.0*tmin) )) + 1
      print*, 'nt =',nt
      do 5 i=1,nt
c        tlog=float(i-1)*dtt+log10(tmin)
        tlog=float(i-1)*dtt+log10(2.0*tmin)
        if (i .eq. 1) then
          tim_wan(i)=10**tlog
c          print*,i,tlog,tim_wan(i)
          tlast=tim_wan(i)
          ntt=1
        else
          t_test=10**tlog
          if (t_test-tlast .ge. tmin) then
            ntt=ntt+1
            nx=int(t_test/tmin)
            t_test=float(nx)*tmin
            tim_wan(ntt)=t_test
            tlast=t_test
c            print*,ntt,i,tlog,tim_wan(ntt)
          end if
        end if
5     continue
      print*,' Number of intervals to compute wander is:', ntt
      if (ntt .gt. max_wan) then
c  test if number of wander periods exceeds dimension; if so, increase spacing
        dtt=dtt*1.5
        go to 4
      end if
c
c  change data spacing 
c
      do 6 i=1,max_flt
6     dat(i)=fmiss
      do 7 i=1,ic
      ix=int((t(i)-t(1))/dt_sam) + 1
      dat(ix)=d(i)
7     continue
      ix_max=int((t(ic)-t(1))/dt_sam) + 1
      if (ix_max .gt. max_flt) then
        print*, 'ix_max=', ix_max,' exceeds max_flt=', max_flt
        stop
      end if


      call wander(tim_wan,wan_dat,max_wan,ntt,dat,max_flt,ix_max,
     &  dt_sam)

c
c  revise list of tim_wan for missing data
c
      nx=ntt
      do 10 i=1,ntt
      if (wan_dat(i) .eq. fmiss) then
        nx=nx-1
        do 11 j=i,nx
          tim_wan(j)=tim_wan(j+1)
          wan_dat(j)=wan_dat(j+1)
11      continue
      end if
10    continue
      if (nx .ne. ntt) then
        print*,' found ',ntt-nx,' periods without estimating wander'
        ntt=nx
      end if
      print*,' Wander**2 of data is'
      do 8 i=1,ntt
8     print*,tim_wan(i),wan_dat(i)
c
c Simulate data using noise model
c
c  read in number of simulations
c
      print*,' Input the number of sets of sythetic data for'
     &  ,' simulations'
      read(5,*)nloop
      write(25,104)nloop
104   format(i5,'     # number of simulations ')
      if (nloop .gt. max_loop) then
         print*,' Number of specified loops exceeds dimensions',
     &    max_loop
         nloop=max_loop
      end if
      print*,' Input the noise model construction'
      print*,'  q -- quadrature add -- convolve each noise filter with'
      print*,'        white noise then add'
      print*,'  n -- same as q'
      print*,'  a -- additive -- add filter functions then convolve'
      print*,'         with white noise'
      read(5,105)ModType
105   format(a1)
      write(25,105)Modtype
c
c  determine filter coeficients for 2 power law noise models, and 1 band-pass filter
c     white noise
c
      print*,' For bandpass filtered noise, input the pass bands;'
      print*,'   flow and fhigh is c/yrs'
      read(5,*)fl,fh
      if (fl .gt. fh) then
        ftemp=fl
        fl=fh
        fh=ftemp
      end if
      print*,'  and input the number of poles'
      read(5,*)np
      write(25,106)fl,fh,np
106   format(2f10.3 '  #  flow fhigh of BP noise c/yr ',/,
     &  i5,'   #  number of poles of BP noise')
c
c  create the filter functions to creat noise
c
      ts=sngl(dt_sam/365.25)
      npts=ix
      do 17 i=1,max_flt
        wh(i)=0.0
        pl1(i)=0.0
        pl2(i)=0.0
        bp(i)=0.0
17    continue      
      if (sig .ne. 0) then
c  white noise impulse response
        wh(1)=1.0d+0
        do 18 i=2,npts+100
          wh(i)=0.0d+0
18     continue
      end if
      if (amp1 .ne. 0) then
c  First power law impulse response

        call frac_diff(pl1,exp1,alpha,ts,npts+100,max_time)
      end if
      if (amp2 .ne. 0) then
c  Second power law impulse response

        call frac_diff(pl2,exp2,0.0,ts,npts+100,max_time)
c        print*,(i,f3(i),i=1,10)
      end if
      if (amp_bp .ne. 0.0) then
        call band_pass_filt(dble(ts),fl,fh,np,max_time,npts+100,bp,1.0)
      end if
      do 30 il=1,nloop
31      continue
        call MakeNoise(dat,Modtype,max_flt,ix,dt_sam,iseed)
c  variance of data
        var=0
        do 38 i=1,npts
38      if (dat(i) .ne. fmiss) var=var+dat(i)**2
        sd=sqrt(var/float(ic))

        print*,'For synthetic data ',il,' Standard deviation is',sd
        iseed=int(dat(1)*100000)
        if (sd .gt. 100000.) go to 31
c  Calculate wander of synthetic data

        call wander(tim_wan,wan_sim,max_wan,ntt,dat,max_flt,ix_max,
     &  dt_sam)
c  save the results
        do 36 i=1,ntt
36      save_wan(il,i)=wan_sim(i)
30    continue

c
c Done with computing wander of synthetic data
c
c  for each period, put the wander in sequential order to do statistics
c
      do 40 j=1,ntt
        do 41 i=1,nloop
41      temp(i)=save_wan(i,j)
        call chron(temp,nloop,max_loop)
        do 42 i=1,nloop
42      save_wan(i,j)=temp(i)
40    continue
      n95=int(nloop*0.025 +0.5)
      n90=int(nloop*0.05 + 0.5)
      n68=int(nloop*0.16 + 0.5)
      n50=int(nloop*0.50 + 0.5)
      print*,' '
      print*,'for wander'
      print*,'68% confidence intervals between ',n68,' and ',nloop-n68
      print*,'90% confidence intervals between ',n90,' and ',nloop-n90
      print*,'95% confidence intervals between ',n95,' and ',nloop-n95
      print*,' median is ',n50
      open(88,file='wander.out')

      RMS_wan=0.
      RMS_wan_long=0.
      RMS_wan_short=0.
      n_wan_short=0
      n_wan_long=0
      n_wan=0.
      do 50 j=1,ntt
        write(88,883)tim_wan(j),sqrt(wan_dat(j)),
     &   sqrt(save_wan(n68,j)),sqrt(save_wan(nloop-n68,j)),
     &   sqrt(save_wan(n90,j)),sqrt(save_wan(nloop-n90,j)),
     &   sqrt(save_wan(n95,j)),sqrt(save_wan(nloop-n95,j)),
     &   sqrt(save_wan(n50,j)),sqrt(wan_var),conf_i*100
c     &   sqrt(wan_mean),sqrt(wan_var),conf_i*100
883     format(f18.9,9f12.4,f8.1)
50    continue
      print*,' '
      print*,' data wander and confidence limits for plotting ',
     &  'in "wander.out" '
      print*,' col 1  is time interval in days, col 2 is ',
     &  'sqrt of data wander'
      print*,' col 3 and 4 is 68% CI, 5 and 6 is 90%, ',
     &   ' and 7 and 8 is 95%'
      print*,'  col 9 is the Median, col 10 is the standard deviation'
      print*,'  col 11 is the confidence that the real data deviates',
     &  ' from simulated data'
      close(88)
      print*,'  Journal file of input in wand_in.jrn'
      close(25)
      rewind(23)
      iseed=int(dat(1)*100000)
      write(23,*)iseed
      close(23)
      stop
      end
      subroutine wander(tim_wan,wan_data,max_wan,ntt,dat,max_flt,ix,
     &  dt_sam)

c  calculates the wander for specificed intervals
c   wan_data(T)=sum (dat(T+t)-dat(t))^2 / N
c  tim_wan are the time interval (days) to estimate wander--input
c  wan_dat are the estimates of wander for specified tim_wan--output
c  dat observations (missing data 99999)
c  max_wan and max_flt  array dimensions
c  ntt  number of periods in tim_wan
c  ic number of data
c  dt_sam sampling interval in days

      dimension tim_wan(max_wan),wan_data(max_wan),
     &  dat(max_flt)
      double precision dt_sam
c  missing data is 99999.
      fmiss=99999.
      do 10 j=1,ntt
        wan_data(j)=0.0
        n=0
        tau=tim_wan(j)
        itau=int(tau/dt_sam)
        if (j .eq. 1) itau_min=int(tim_wan(1)/dt_sam)
        if (j .gt. 1) itau_min=int((tim_wan(j-1))/dt_sam +1)
        itau_max=int(tim_wan(j)/dt_sam)
c        print*,tim_wan(j),itau_min,itau_max
        do 151 k_tau=itau_min,itau_max
          do 152 i=1,ix-k_tau
            if ((dat(i) .ne. fmiss ) .and. 
     &        (dat(i+k_tau) .ne. fmiss)) then
              n=n+1
             wan_data(j)=wan_data(j)+(dat(i+k_tau)-dat(i))**2
           end if
152       continue
151     continue
        if (n .ne. 0) wan_data(j)=(wan_data(j)/float(n))
        if (n .eq. 0) wan_data(j)=fmiss
c      print*,j,tim_wan(j),itau,itau_min,itau_max,n,wan_data(j)

10    continue
      return
      end
      subroutine MakeNoise(ran,Modtype,max_time,npts,dt_sam,iseed)
      dimension ran(max_time),ran1(131072),ran2(131072),ran3(131072),
     &  ran4(131072)
      double precision dt_sam,f1(131072),f2(131072),
     & f3(131072),f4(131072),rnx,ZBQLNOR,filt(131072)
      character*1 Modtype
      common /NoiseMod/expon1,amp1,alpha,expon2,amp2,sigma,
     & fl,fh,np,f_amp,fmiss
      common /Filter/f1,f2,f3,f4
c      print*,expon1,amp1,alpha,expon2,amp2,sigma,
c     & fl,fh,np,f_amp,fmiss

      ts=sngl(dt_sam/365.25)
      dts=(dt_sam/365.25)
c      print*,Modtype,max_time,npts,dt_sam,ts
      do 9 i=1,npts+100
        ran1(i)=0.0
        ran2(i)=0.0
        ran3(i)=0.0
        ran4(i)=0.0

9     continue

c      print*,'iseed',iseed

      call srand(iseed)
      if ((ModType .eq. 'n') .or. (ModType .eq. 'q')) then
c        print*,' standard noise of 4 independent noise models'
        if (sigma .ne. 0) then
          do 10 i=1,npts+100
          ran1(i)=sigma*sngl(ZBQLNOR(0.0d+0,1.0d+0))
10        continue
          call convolv(npts+100,f1,ran1,max_time)
        end if
        if (amp1 .ne. 0) then
          do 11 i=1,npts+100
          ran2(i)=amp1*sngl(ZBQLNOR(0.0d+0,1.0d+0))*((ts)**(expon1/4))
11        continue
          call convolv(npts+100,f2,ran2,max_time)
        end if
        if (amp2 .ne. 0) then
          do 12 i=1,npts+100
          ran3(i)=amp2*sngl(ZBQLNOR(0.0d+0,1.0d+0))*((ts)**(expon2/4))
12        continue
          call convolv(npts+100,f3,ran3,max_time)
        end if

        if (f_amp .ne. 0) then
          do 13 i=1,npts+100
          ran4(i)=f_amp*sngl(ZBQLNOR(0.0d+0,1.0d+0))
13        continue
          call convolv(npts+100,f4,ran4,max_time)
        end if
        do 15 i=101,npts
          ran1(i)=ran1(i)+ran2(i)+ran3(i)+ran4(i)
15      continue
      else
c        print*,' sum impulse responses and then add'
        do 20 i=1,npts+100
          filt(i)=f1(i)*sigma+amp1*f2(i)*((ts)**(expon1/4))+
     &       amp2*f3(i)*((ts)**(expon2/4))+f_amp*f4(i)

          ran1(i)=1.0*sngl(ZBQLNOR(0.0d+0,1.0d+0))
20      continue 


        call convolv(npts+100,filt,ran1,max_time)

      end if
c
c   Resample the noise at time of measurement
c

        do 58 i=1,npts
          if (ran(i) .ne. fmiss) ran(i)=ran1(i+100)


58      continue


      return
      end
      subroutine convolv(npts,H,ran,max)
      dimension ran(max),ranin(max)
      double precision H(max)
      do 40 i=1,npts
      ranin(i)=ran(i)
40    continue
      do 52 i=1,npts
        ran(i)=0.0
        do 53 j=1,i
        ran(i)=ran(i)+ranin(j)*H(i-j+1)
53      continue
52    continue
      return
      end
      subroutine chron(t,ic,max)
c  puts data into chronological order
c  code riped-off from numerical recipes piksrt.f
      real t(max),tx
      do 12 j=2,ic
         tx=t(j)

         do 11 i=j-1,1,-1
         if (t(i) .eq. tx) t(i)=t(i)+0.001
         if (t(i) .lt. tx) go to 10
           t(i+1)=t(i)
11       continue
         i=0
10       t(i+1)=tx
12    continue
      return
      end
