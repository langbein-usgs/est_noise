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
