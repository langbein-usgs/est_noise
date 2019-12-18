c  Program computes PSD curve for given input of PSD functions
c    Functions include
c      white noise
c      power law noise
c      power law noise with Gauss-Markov roll-off
c      band pass filtered
c      spiked noise
c  MODIFICATION of psd_cal6 to include computation of PSD
c   using quadrature addition (which is psd_calc) and addition
c   of filter function (as est_noise7.2x using modtype a, f, or c)
c    aka simple addition
c
      double precision ran1(131072), ran2(131072),ran3(131072),
     & ran4(131072),ran(131072)
      double precision wsave(393231),dran(131072),
     & dreal(131072),ddimag(131072),azero
c  swave dim 3*n+15
      character*1 modtype
      nmax=131072
      max=131072
      pi=3.14159265
      open(10,file='psdcalc.out')
      print*,' Input the sampling interval in days'
      read(5,*)tsam
      fny=0.5*(365.25/tsam)
      fs=365.25/tsam
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
      fgm=fgm/tsam
c      psd_db=2.8384+20.*alog10(amp_pl) - 8.0*fexp
c      psd_db=5.0*fexp1*alog10(fs)-7.399*fexp1-10*alog10(fs)+2.076
      psd_db=10.0*alog10(2.0)-10.*fexp1*alog10(2.0*3.1415926)
     & +5.0*fexp1*alog10(fs)-10.*alog10(fs)
c      print*,psd_db,psd_dbx
      
      if (abs(amp_pl) .ne. 0) then
         pow_db=20*alog10(amp_pl)+psd_db
         print*,' first power law',pow_db,' db', amp_pl,fexp 
         P1o=10**(pow_db/10.)
      else
         P1o=0.0
      end if
c      print*,psd
      print*,' Input parameters for SECOND power law noise'
      print*,' Input the index (n=2 for random walk noise)'
      read(5,*)fexp2
      print*,' Input the amplitude'
      read(5,*)amp_pl
      amp_pl2=amp_pl
c      if (amp_pl .eq. 0) amp_pl=0.00001
c      psd_db=2.8384+20.*alog10(amp_pl) - 8.0*fexp
     
c      psd_db=5.0*fexp2*alog10(fs)-7.399*fexp2-10*alog10(fs)+2.076
      psd_db=10.0*alog10(2.0)-10.*fexp2*alog10(2.0*3.1415926)
     & +5.0*fexp2*alog10(fs)-10.*alog10(fs)      
c      print*,5.0*fexp2*alog10(fs),7.399*fexp2,10*alog10(fs),2.076
c      psd_db1=5.0*fexp2*log10(fs)-7.399*fexp2-10*log10(fs)+2.076
c      print*,fexp2,fs,psd_db,pow_db1
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
      print*,Pw,10.*alog10(Pw)
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
      print*,' '
c      print*,' Input spike noise'
c      print*,' Input its frequency, c/yr'
c      read(5,*)fspike
c      print*,' Input its amplitude'
c      read(5,*)amp_sp
c      Ps=(1.0/fl)*(amp_sp**2)/2.

c      print*,P1o,fexp1,P2o,fexp2
c  
c  convert from radians/yr to cycles/yr
       fgm=fgm/(2.0*3.14159265)
      print*," "
      print*,"  Input the noise type using n, q, or a"
      print*,"    n or q -- quadrature addition"
      print*,"    a -- simple addition"
      read(5,105)modtype
105   format(a1)
      if ((modtype .eq. "n" ) .or. (modtype .eq. "q")) then
      print*,"  PSD using quadrature"
      do 1 i=1,npts
        freq=fl*float(i)
c  power law
        powPL1=P1o
         powPL2=P2o

c  modification per simon williams
       if (fexp1 .ne. 0) powPL1=P1o/((fgm**2 + freq**2)**(fexp1/2.0))
c       if (fexp1 .ne. 0) powPL1=P1o/(fgm**fexp1 + freq**fexp1)

       if (fexp2 .ne. 0) powPL2=P2o/( freq**fexp2)
c        powPL1=P1o/(fgm**fexp1 + freq**fexp1)
c        powPL2=P2o/( freq**fexp2)

c  white noise 
        Pw=Pw  

c bandpass filtered
        h=((freq/flow)**2)
     &    /((1.0+((freq/flow)**2))*(1.0+((freq/fhigh)**2)))
c        h=(amp_bp**2)*h
        if (npole .eq. 2) h=h*h
        if (npole .eq. 3) h=h*h*h
        if (npole .eq. 4) h=h*h*h*h
        if (npole .eq. 5) h=h*h*h*h*h
        if (npole .eq. 6) h=h*h*h*h*h*h

        powBP=(4.0*amp_bp**2)*h
c        print*,freq,h

c spike noise
c        df=abs(freq-fspike)
        powSP=0.
c        if (df .lt. 0.5*fl) powSP=Ps
       
c Sum the noise
        
        pow=powPL1+Pw+powBP+powSP+powPL2
        power=10.*alog10(pow)

c      write(11,*)freq,pow,powPL1,Pw,powBP,powSP,powPL2
c      write(9,*)freq,power, 10.*alog10(powPL1),10.*alog10(Pw),
c     & 10.*alog10(powPL2)       
      write(10,*)freq,power
1     continue 
      else
        print*,"  compute PSD for simple addition using fft"
c  Computes filter function for each of the above, sums the each filter function
c   then computes the FFT of the filter function and converts the FFT to PSD by dividing by dfreq.

      fl=1/tlen
      npts=int(fs/fl)+1
c
c  force npts be a power of 2
c
      nn=int(alog(float(npts))/alog(2.))+1
c      print*,' nn=',nn
c  double the required amount since time series gets shifted by original time-series length and padded with 0.
      npts=2*(2**nn)
      if (npts .gt. nmax) then
        print*,' Length of time-series with padding exceeds dimensions'
        print*,' npts=', npts,' nmax=',nmax
        stop
      end if
      fl=(365.25/tsam)/float(npts)
      tsam=tsam/365.25
c      print*," tsam= ",tsam," fl= ", fl
      if (sigma .ne. 0.) then
        call white(nmax,sigma,ran1,npts)
      end if
      if (amp_pl1 .ne. 0.0) then
c        alpha=fgm/(2.0*pi)
c  convert from cycles/yr to radians/yr 
        alpha=fgm*(2.0*pi)
c        print*,nmax,npts,iseed,tsam,fexp1,amp_pl1,alpha
c        call power_law(nmax,npts,iseed,tsam,ran2,fexp1,amp_pl1,alpha)
        call frac_diff(ran2,fexp1,alpha,tsam,npts,nmax)
        amp_pl1=amp_pl1*(tsam**(fexp1/4))

        do 20 i=1,npts
        ran2(i)=ran2(i)*amp_pl1
20      continue
      end if
      if (amp_pl2 .ne. 0.0) then
        alpha=0.0
c        call power_law(nmax,npts,iseed,tsam,ran3,fexp2,amp_pl2,alpha)
        call frac_diff(ran3,fexp2,alpha,tsam,npts,nmax)
        amp_pl2=amp_pl2*(tsam**(fexp2/4))
        do 21 i=1,npts
        ran3(i)=ran3(i)*amp_pl2
21      continue
      end if
      if (amp_bp .ne. 0.) then

        call band_pass_filt(dble(tsam),flow,fhigh,npole,nmax,npts,
     & ran4,amp_bp)
      end if
      do 19 i=1,nmax
      ran(i)=0.0
19    continue
c Add up the sequencies of random noise

      if (sigma .ne. 0.0) call sum_ran(nmax,npts,ran1,ran)
      if (amp_pl1 .ne. 0.0) call sum_ran(nmax,npts,ran2,ran)
      if (amp_pl2 .ne. 0.0) call sum_ran(nmax,npts,ran3,ran)
      if (amp_bp .ne. 0.0) call sum_ran(nmax,npts,ran4,ran)

c
c  shift time series by npts points
c

      ishift=npts/2
      nx=npts-ishift
c      print*,' ishift=',ishift,nx

      do 28 i=nx,1,-1
c        write(71,*)i,ishift+i,ran(i)
        ran(ishift+i)=ran(i)
c        write(71,*)ishift+i,ran(ishift+i)
28    continue
      do 281 i=1,ishift
        ran(i)=0.0d+0
        dran(i)=0.0d+0
281   continue
c      print*,(ran(i),i=ishift,ishift+10)
c
c  apply a cosine taper to last portion of filter
      ntap=ishift/1
c      do 282 i=1,ishift
      do 282 i=1,ntap
        taper=cos(2*pi*float(i)/(4.0*float(ntap)))

c        print*,i,npts-ishift+i,ran(npts-ishift+i),taper
        dran(npts-ntap+i)=dble(taper*ran(npts-ntap+i))
282   continue
      print*,' Output the filter function  filt.dat'
      open(7,file='filt.dat')

      do 30 i=1,npts
        write(7,*)dran(i)
c        write(71,*) ran(i)
30    continue
      close (7)
c
c  Compute FFT
c
c  initialize
c      print*,' initialize'
      call dzffti(npts,wsave)
c  do FFT
      call dzfftf(npts,dran,azero,dreal,ddimag,wsave)
      do 29 i=2,(npts/2)-1
c        ran(i)=ran(i)/float(npts/2)
c        ran(i)=ran(i)
c        write(81,*)i,dreal(i),ddimag(i)
        Pdb=10.0*dlog10((dreal(i)**2 + ddimag(i)**2)/fl)
     &        + 10.0*alog10(float(npts)/2.)
          write(10,*)fl*float(i),Pdb
29    continue
      i=(npts/2)
c      write(81,*)i,dreal(i),ddimag(i)
      Pdb=10.0*dlog10((dreal(i-1)**2 + ddimag(i)**2)/fl)
     &        + 10.0*alog10(float(npts)/2.)
      write(10,*)fl*float(i),Pdb
      end if  
      print*,'  Computed PSD in psdcalc.out'
      close(10)     
      stop
      end
      subroutine sum_ran(nmax,npts,ran1,ran)
c
c Adds time series
c  
c  nmax is the dimension of data array (input)
c  npts  is number of points (input)
c  ran1 is input time series
c  ran is replaced by ran+ran1
      double precision ran(nmax),ran1(nmax)
      do 1 i=1,npts
1     ran(i)=ran(i)+ran1(i)
      return
      end     
      subroutine white(nmax,sigma,ran,npts)
c  creates a filter function for white noise with 
c  nmax is the dimension of data array (input)
c  sig  is the standard deviation of white noise (input)
c  ran  filter function for white noise (output)
c  npts  is length of time series number of points (input)
c  iseed is the seed for random number generator (input/output)
      double precision ran(nmax)

      do 1 i=1,npts
1     ran(i)=0.0d+0
      ran(1)=sigma
      do 2 i=npts,nmax
2     ran(i)=0.0d+0
 
      return
      end      
