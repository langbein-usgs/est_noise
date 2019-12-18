c  applies offsets, rate changes, and exponential to GPS data
c
c  Modification of existing program off_rate_exp_gps_2
c    to accommodate different formats of time stamps
      dimension yr(3009901), fjul(3009901), data(3009901),
     &   er(3009901),time(3009901)
      double precision time,toff,ffjul,ffjul2, toff2,fjul,timed,tref
     &  ,dec_timed
      character*16 ifile, ofile
      character*1 func_type,sgn
      character*3 netd, nett
      character*80 string
      max=3009901
      print*,' Program performs adjustments to GPS data'
      print*,'   works with offset, rate changes and exponentials'
      print*,'  specify function type'
      print*,'   followed by specifications for function'
      print*,'   Specification of each function covers several lines'
      print*,'  o = offset'
      print*,'     time'
      print*,'     ammount of offset'
      print*,'  R = rate'
      print*,'      rate in unit/yr'
      print*,'  r = rate change'
      print*,'     time1, time2'
      print*,'     rate change unit/yr'
      print*,'  e = exponential (1 - e (-t/tau))'
      print*,'     time'
      print*,'      tau (years), amplitude'
      print*,'  m = logrithm  log10 (1 + t/tau) '
      print*,'     time'
      print*,'      tau (years), amplitude'
      print*,'  s = cosine/sine '
      print*,'     reference time'
      print*,'  , period (days)  cosine amp, sine amp'
      print*," "

      print*,'  Input the data format type'
      print*,' otr=data'
      print*,' otr   format of year day_of_year obs err'

      print*,' otd   format of YearMnDa  obs err'

      print*,' otx   format of year mo da obs err'

      print*,' mjd   Modified Julian day, obs, err'
      print*,'For all of the above, the da and doy are ',
     &  'double precision'
      print*,'    which allows decimal days to accomodate seconds'
      print*,' '
      print*,' gmt  GMT format for time'
      print*,'        year-mo-daThr:mn:secs.x  obs  err'
      read(5,105)netd
105   format(a3)
      print*,' Input GPS data file to apply offsets and ',
     & 'adjustments (cmbb.u)'
      read(5,100)ifile
100   format(a16)
      open(1,file=ifile)
      print*,' Input the format of times used in adjustment file'
      print*,'    See above listing'
      read(5,105)nett
      print*,' Input name of file with offset corrections (Ocmbb.u)'
      read(5,100) ofile
      print*,'  file name is  ',ofile
      print*,'  Input either d/a to indicate whether adjustments'
      print*,'   will (d)etrend data, or (a)dd a trend to data'
      print*,'   a d will multiply given adjustments by -1.0'
      read(5,102)sgn
102   format(a1)
      sn=1.0
      if (sgn .eq. 'd') sn=-1.0

      open(2,file=ofile)
      n=1
1     continue
        if (n .gt. max ) then
          print*,' Number of data',n,' exceeds dimension',max
          stop
        end if
        read(1,106,end=3) string
106     format(a80)
        call GetData(netd,string,timed,obs,err)
        data(n)=obs
        er(n)=err
        time(n)=timed
  
        n=n+1
        go to 1
3     continue
      n=n-1
      if (n .gt. max) then
        print*,' number of data exceeds dimensions', n,max
        stop
      end if
      print*,' number of data',n
      noff=0
      nexp=0
      nrate=0
      nlog=0
      nsea=0
      nrates=0
20    continue
      read(2,*,end=50) func_type
        if (func_type .eq. 'o') then
c    deal with offset
          call GetTime3(nett,toff)
          read(2,*)  off
c          print*,' Offset at ',toff
          do 21 i=1,n
          if (time(i) .ge. toff) data(i)=data(i)+sn*off
21        continue
          noff=noff+1
        end if
        if (func_type .eq. 'r') then
c  deal with rate change
          call GetTime2(nett,toff,toff2)
c          print*,' rate change between', toff,' and ',toff2
          read(2,*) rate
          do 22 i=1,n
          if ((time(i) .ge. toff) .and. 
     &     (time(i) .lt. toff2))  data(i)=data(i)
     &      +sn*rate*(time(i)-toff)/365.25
          if (time(i) .ge. toff2) data(i)=data(i)
     &      + sn*rate*(toff2-toff)/365.25
22        continue

          nrate=nrate+1
        end if
        if (func_type .eq. 'e') then
c  deal with exponentials
          call GetTime3(nett,toff)
c          print*,' exponential at ',toff
          read(2,*) tau, amp
          do 23 i=1,n
          if (time(i) .ge. toff) data(i)=data(i)
     &     +sn*amp*(1.0 - exp(-(time(i)-toff)/(tau*365.25)) )
23        continue
          nexp=nexp+1
        end if
        if (func_type .eq. 'm') then
c  deal with logrithm
          call GetTime3(nett,toff)
 
          read(2,*) tau, amp
c         print*,' log function at ',toff
          do 24 i=1,n
          if (time(i) .ge. toff) data(i)=data(i)
     &     +sn*amp*dlog10(1.0 + (time(i)-toff)/(tau*365.25))
24        continue
          nlog=nlog+1
        end if
        if (func_type .eq. 's') then
          call GetTime3(nett,tref)

          read(2,*) per,camp,samp
          freq=2.0*3.1415926/per
          do 25 i=1,n
          data(i)=data(i)+sn*camp*cos(freq*(time(i)-tref))
     &      + sn*samp*sin(freq*(time(i)-tref))
25        continue
         nsea=nsea+1
        end if 
        if (func_type .eq. 'R' ) then
c  deal with removing the rate
          read(2,*)rate
          do 26 i=1,n
          data(i)=data(i)+sn*rate*(time(i)-time(1))/365.25
26        continue
          nrates=nrates+1
       end if
      go to 20
50    continue
      rewind(1)
      do 60 i=1,n
        itime=int(time(i))
        timed=(time(i))-int(time(i))
        dec_timed=timed
        if (netd .eq. 'otr' ) then
          call inv_jul_time(itime,nyr,jul)
          write(1,6501)nyr,float(jul)+dec_timed,data(i),er(i)
        end if
        if (netd .eq. 'otx' ) then
          call invcal(itime,nyr,mn,idate)
          write(1,6502)nyr,mn,idate+dec_timed,data(i),er(i)
        end if
        if (netd .eq. 'otd' ) then
          call invcal(itime,nyr,mn,idate)
          write(1,6503)nyr,mn,idate,dec_timed,data(i),er(i)
        end if 
        if (netd .eq. 'gmt' ) then
          call invcal(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_timed)
          imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
          sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          write(1,6504)nyr,mn,idate,ihr1,imn1,isec1,sec1,data(i),er(i)
        end if
        if (netd .eq. 'mjd' ) then
          write(1,6505)dble(itime)+dec_timed+36933.d+0,data(i),er(i)
        end if
60    continue
6501  format(1x,i5,1x, f15.9, 1x, f15.2,1x,f8.2)
6502  format(1x,2i5,f15.9, 1x,f15.2,1x,f8.2 )
6503  format(1x,i4,i2.2,i2.2,f10.9, 1x,f15.2,1x,f8.2)
6504  format(1x,i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1, 1x,f15.2,1x,f8.2)
6505  format(1x,f18.9, 1x,f15.2,1x,f8.2)
      print*,' number of sinusoids fixed', nsea
      print*,' number of offset fixed', noff
      print*,' number of ratechanges fixed', nrate
      print*,' Number of exponentials fixed', nexp
      print*,' Number of secular rates fixed', nrates
      print*,' Number of logrithm/omori fixed', nlog
      print*,' Adjusted data in file ', ifile
      close(1)
      close(2)
      stop
      end
      subroutine GetTime2(net,t_start,t_stop)
c
c   read time interval to estimate a RATE change
c
c  net is input for the time format
c  t_start and _stop are output: time span, in days since 1960, to analyze data
c
      character*3 net
      character*80 string1,string2
      double precision fjul1, fjul2, t_start, t_stop,day1,day2
      if (net .eq. 'otr') then
c        print*,'rate change interal',n,' year day_of_yr year day_of_yr'
c        print*,'  "day_of_yr" may be decimal day'
        read(2,*)yr1,fjul1,yr2,fjul2
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        t_start=dble(iday)+fjul1-int(fjul1)

        call jul2day(int(yr2),int(fjul2),iday)
c days since 1960
        t_stop=dble(iday)+fjul2-int(fjul2)
      end if 
      if (net .eq. 'otx' ) then
c        print*,'  rate change interal',n,' year mn da  year mn da'
c        print*,'   "da" may be a decimal day'
        read(2,*) yr1,mn1,day1,yr2,mn2,day2
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        call date2day(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'otd' ) then
c        print*,' rate change interal',n,' YearMmnDa  YearMmnDay'
c        print*,'   "da" may be a decimal day'
        read(2,*) day1, day2
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
c        print*,'  rate change interal',n,' MJD  MJD'
c        print*,'   May be decimal days'
        read(2,*) day1, day2
        t_start=day1-36933.0
        t_stop=day2-36933.0
      end if
      if (net .eq. 'gmt' ) then
c        print*,'  rate change interal',n,' in GMT time format'
c        print*,'   two entries required'
        read(2,*) string1, string2
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
      subroutine GetTime3(net,t_start)
c
c   read time to estimate an OFFSETs
c
c  net is input for the time format
c
c  t_start day since 1960

c
      character*3 net
      character*80 string1
      double precision fjul1, t_start,day1
      if (net .eq. 'otr') then
        read(2,*)yr1,fjul1
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        t_start=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        read(2,*) yr1,mn1,day1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        read(2,*) day1
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        read(2,*) day1
        t_start=day1-36933.0
      end if
      if (net .eq. 'gmt' ) then

        read(2,*) string1
        read(string1,100)iyr,imn,ida,ihr,mn,sec
100     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)
        call date2day(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)
      end if
      return
      end


