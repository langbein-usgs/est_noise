c  generates a time series of correlated noise.
       double precision dt_sam, t(1531072),dec_timed
       dimension d(1531072)
       character*1 ModType
       character*3 net
       max_time=1531072
      open(1,file="seed.dat")
       read(1,*)iseed
      open(25,file='gen_in.jrn')

       print*,'  Program generates power law noise; more to follow'
       print*,' Input the method to create power law noise'
       print*,'  n --  noise added from independent sources ',
     &  '(quadrature)'
       print*,'  a -- filter functions added than convolved',
     &  '    with single noise source "additive"'
       read(5,100)ModType
100    format(a1)
       write(25,100)Modtype
       print*,' Input the sampling interval in days'
       print*,'   1= 1 day'
       print*,'   0.0416667= 1 hour'
       print*,'   0.0069444444= 10 minutes'
       print*,'   0.00069444444= 1 minute'
       print*,'   1.157407e-05 = 1 second'
       read(5,*)dt_sam
       write(25,101)dt_sam
101    format(e16.7,'     # sampling interval in days')
       print*,' Length of time series in years'
       read(5,*)t_len
       write(25,102)t_len
102    format(f16.7,'     # time series length in years')
       npts=int(t_len*365.25/dt_sam)
       print*,' Number of point in time-series'
       if (npts .gt. max_time) then
         print*,' Number of points', npts,' exceeds dimension',max_time
         stop
       end if
       do 1 i=1,npts
         t(i)=(float(i-1))*dt_sam
         d(i)=0.0
1      continue
       print*,' What time format should the noisy data have?'
       print*,'   otr --  year, day_of_year'
       print*,'   otx --  year Mn Da'
       print*,'   otd --  yearMnDay'
       print*,'   gmt --  GMT compatible time'
       print*,'   mjd --  modified julian day'
       read(5,103)net
103    format(a3)
       write(25,103)net
       call genNoise(ModType,t,d,npts,max_data,iseed,
     &  dt_sam,max_time)
       rewind(1)
       write(1,*)iseed
       close(1)
       close(25)
       print*,' Journal file of input in gen_in.jrn'
       open(2,file='gen.out')
       open(3,file='noise.dat')
       do 10 i=1,npts
         t(i)=t(i)+1.0
         write(2,*)d(i)
         itime=int(t(i))
         dec_timed=(t(i))-int(t(i))
         if (net .eq. 'otr' ) then
           call inv_jul_time(itime,nyr,jul)
           write(3,6501)float(nyr),float(jul)+dec_timed,d(i)
         end if
         if (net .eq. 'otx' ) then
           call invcal(itime,nyr,mn,idate)
           write(3,6502)nyr,mn,idate+dec_timed,d(i)
         end if
         if (net .eq. 'otd' ) then
           call invcal(itime,nyr,mn,idate)
           write(3,6503)nyr,mn,idate,dec_timed,d(i)
         end if 
         if (net .eq. 'gmt' ) then
           call invcal(itime,nyr,mn,idate)
c           write(3,*)i,t(i),itime,nyr,mn,idate
           ihr1=int(24.0*dec_timed)
           imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
           sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
           sec1=sec1*3600.0*24.0
           isec1=int(sec1)
           write(3,6504)nyr,mn,idate,ihr1,imn1,isec1,sec1,d(i)
         end if
         if (net .eq. 'mjd' ) then
            write(3,6505)dble(itime)+dec_timed+36933.d+0,d(i)
         end if       
10     continue
6501  format(1x,f6.0,1x, f15.9, 1x, f15.2, ' 1.0')
6502  format(1x,2i5,f15.9, 1x,f15.2, '  1.0')
6503  format(1x,i4,i2.2,i2.2,f10.9, 1x,f15.2, ' 1.0')
6504  format(1x,i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1, 1x,f15.2, ' 1.0')
6505  format(1x,f18.9, 1x,f15.2, ' 1.0')
       close(3)
       close(2)
       print*,' Sequence of correlated data in gen.out; with ',
     &  'no time stamps'
       print*,' Simulated noisy data in noise.dat'
       stop
       end
