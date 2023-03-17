      subroutine genNoise(ModType,t,d,npts,max_data,iseed,
     &  dt_sam,max_time)
c
c  substitutes noise for real data
c
c  ModType, n, a, c, f
c         n and q is standard, noise summed from independent random number generators
c         a, c, or f sums the transfer function and uses single random number generator
c  t is time
c  d is data-  which will be filled with noise
c  ic and max_data are number of data and the dimensions of arrays
      dimension d(max_data),ran1(32768),ran2(32768),ran3(32768),
     &  ran4(32768)
      double precision t(max_data),dt_sam,f1(32768),f2(32768),
     & f3(32768),f4(32768),rnx,ZBQLNOR
c      dimension d(max_data),ran1(2097151),ran2(2097151),ran3(2097151),
c     &  ran4(2097151)
c      double precision t(max_data),dt_sam,f1(2097151),f2(2097151),
c     & f3(2097151),f4(2097151),rnx,ZBQLNOR
      character*1 ModType
 
        print*,' Program creates colored noise with'
        print*,'  PSD   Po / ((alpha^n + f ^ n)) Plus white noise'
        print*,'     plus more power-law noise of P1/f^j'
        print*,'     pluse band-passed filterd noise'
        print*,' '
        print*,' First power law/Gauss-Markov parameters'
        print*,' input the exponent n'
        read(5,*)expon1
        write(25,301)expon1
301     format(f15.3,'     # input noise 1 exp term')
        print*,' input term proportional to Po'
        read(5,*)amp1
        write(25,302)amp1
302     format(f15.3,'     # input noise 1 amp. term')
        print*,' input the time constant alpha in radians/year'
        read(5,*)alpha
        write(25,303)alpha
303     format(f15.3,'     # input noise 1 GM term')
        print*,' '
        print*,' Second power law parameters'
        print*,' input the exponent n'
        read(5,*)expon2
        write(25,304)expon2
304     format(f15.3,'     # input noise 2 exp term')
        print*,' input term proportional to Po'
        read(5,*)amp2
        write(25,305)amp2
305     format(f15.3,'     # input noise 2 amp. term')
        print*,' '
        print*,' Input the white noise'
        read(5,*)sig
        write(25,306)sig
306     format(f15.3,'     # input white noise')
        print*,' '
        print*,' Creates a second set of noise consisting of ',
     & 'band-passed'
        print*,'  band passed filtered, white noise'
        print*,' Input the lower and upper limits of the filter',
     &' in c/yr'
        read(5,*)fl,fh
        if (fl .gt. fh) then
           f_tem=fl
           fl=fh
           fh=f_tem
        end if
        write(25,*)fl,fh
        print*,' Input the number of pole of the filter, 1,2,..'
        read(5,*)np
        write(25,*)np
        print*,' Input the amplitude of the noise'
        read(5,*)f_amp
        write(25,307)f_amp
307     format(f15.3,'     # Input Band Pass amp')     
      ispinUp=100
      ts=sngl(dt_sam/365.25)
      dts=(dt_sam/365.25)
      do 9 i=1,npts+ispinUp
        ran1(i)=0.0
        ran2(i)=0.0
        ran3(i)=0.0
        ran4(i)=0.0
        f1(i)=0.0d+0
        f2(i)=0.0d+0
        f3(i)=0.0d+0
        f4(i)=0.0d+0
9     continue

      if (sig .ne. 0) then
c  white noise impulse response
        f1(1)=1.0d+0
        do 8 i=2,npts+ispinUp
          f1(i)=0.0d+0
8      continue
      end if
      if (amp1 .ne. 0) then
c  First power law impulse response

        call frac_diff(f2,expon1,alpha,ts,npts+ispinUp,max_time)
      end if
      if (amp2 .ne. 0) then
c  Second power law impulse response

        call frac_diff(f3,expon2,0.0,ts,npts+ispinUp,max_time)
c        print*,(i,f3(i),i=1,10)
      end if
      if (f_amp .ne. 0.0) then
        call band_pass_filt(dble(ts),fl,fh,np,max_time,npts+ispinUp,
     &    f4,1.0)
      end if
c
c  Depending upon noise model type, either add all of the implus responses
c    and convolve with noise, or convolve noise with each impulse response, and then add
c
      print*,' ModType ',ModType

      call srand(iseed)
      if (ModType .eq. 'n' ) then
        print*,' standard noise of 4 independent noise models'
        if (sig .ne. 0) then
          do 10 i=1,npts+ispinUp
          ran1(i)=sig*sngl(ZBQLNOR(0.0d+0,1.0d+0))
10        continue
          call convolv(npts+ispinUp,f1,ran1,max_time)
        end if
        if (amp1 .ne. 0) then
          do 11 i=1,npts+ispinUp
          ran2(i)=amp1*sngl(ZBQLNOR(0.0d+0,1.0d+0))*((ts)**(expon1/4))
11        continue
          call convolv(npts+ispinUp,f2,ran2,max_time)
        end if
        if (amp2 .ne. 0) then
          do 12 i=1,npts+ispinUp
          ran3(i)=amp2*sngl(ZBQLNOR(0.0d+0,1.0d+0))*((ts)**(expon2/4))
12        continue
          call convolv(npts+ispinUp,f3,ran3,max_time)
        end if
        if (f_amp .ne. 0) then
          do 13 i=1,npts+ispinUp
          ran4(i)=f_amp*sngl(ZBQLNOR(0.0d+0,1.0d+0))
13        continue
          call convolv(npts+ispinUp,f4,ran4,max_time)
        end if
        do 15 i=1,npts+ispinUp
          ran1(i)=ran1(i)+ran2(i)+ran3(i)+ran4(i)
15      continue
      else
        print*,' sum impulse responses and then add'

        do 20 i=1,npts+ispinUp
          f1(i)=f1(i)*sig+amp1*f2(i)*((ts)**(expon1/4))+
     &       amp2*f3(i)*((ts)**(expon2/4))+f_amp*f4(i)

          ran1(i)=1.0*sngl(ZBQLNOR(0.0d+0,1.0d+0))
20      continue   

        call convolv(npts+ispinUp,f1,ran1,max_time)
c        do 223 i=1,10
c223     print*,i,ran1(i+100)
      end if
c
c   Resample the noise at time of measurementsc
c

        do 58 i=1,npts
          ix=int(t(i)-(t(1))/dt_sam)+1
          d(i)=ran1(i+ispinUp)

58      continue

      iseed=int(1000000.*rand())
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
