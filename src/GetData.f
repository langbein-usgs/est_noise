      subroutine GetData(net,string,time,obs,err)
c
c   Parses data from string depending on net variable
c  outputs timex, observation, and its assigned error (which
c    isn't used in the code, at least not yet)
      character*132 string
      double precision time,fjul1,day1
      character*3 net
      character*80 string1,string2,string3
c      print*,"net= ",net
c      print*,string
      if (net .eq. 'otr') then

        read(string,*)yr1,fjul1,obs,err
c        print*,yr1,fjul1,obs,err
        call jul2day(int(yr1),int(fjul1),iday)
c days since 1960
        time=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        read(string,*)yr1,mn1,day1,obs,err
        call date2day(int(yr1),int(mn1),int(day1),iday)
        time=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        read(string,*)day1,obs,err
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2day(int(yr1),int(mn1),int(day1),iday)
        time=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        read(string,*) day1,obs,err
        time=day1-36933.0
      end if
      if (net .eq. 'gmt' ) then
        read(string,*) string1,string2,string3
        read(string1,100)iyr,imn,ida,ihr,mn,sec
100     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)
        call date2day(iyr,imn,ida,iday)
        time=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) +
     &    dble(sec)/(24.0*3600.0)
        read(string2,*)obs
        read(string3,*)err
      end if
      return
      end
