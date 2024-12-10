module GetStuff_mod

!  Data and time format subroutines

  use time_mod

  implicit none
  private
 
  public ::  GetData, GetTime1, GetTime2, GetTime3, GetTime4, DataOut
  
  contains 
    subroutine GetData(net,string,time,obs,err)
!
!  Parses data from string depending on net variable
!  outputs timex, observation, and its assigned error (which
!    isn't used in the code, at least not yet)  
      use iso_fortran_env
      character(len=132) string
      real(kind=real64)  :: time,fjul1,day1
      character(len=3) :: net
      character(len=80) :: string1,string2,string3
      real(kind=real32) :: yr1,err,obs,sec
      integer :: mn1,iday,mn,iyr,imn,ida,ihr
      if ((net .eq. 'doy') .or. (net .eq. 'otr')) then

        read(string,*)yr1,fjul1,obs,err
        call doy2num(int(yr1),int(fjul1),iday)
                    
! days since 1800
        time=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        read(string,*)yr1,mn1,day1,obs,err
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        time=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        read(string,*)day1,obs,err
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        time=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        read(string,*) day1,obs,err
        time=day1
      end if
      if (net .eq. 'gmt' ) then
        read(string,*) string1,string2,string3
        read(string1,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec
        call date2numJL(iyr,imn,ida,iday)
        time=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)
        read(string2,*)obs
        read(string3,*)err
      end if    
     end subroutine GetData
     
      subroutine GetTime1(net,t_start,t_stop,unitnum,wout)
!
!   read time interval to analyze data
!
!  net is input for the time format
!  t_start and _stop are output: time span, in days since 1800, to analyze data
!  wout --- output date on unit 25
!
      use iso_fortran_env
      character(len=3) :: net
      character(len=80) :: string1,string2
      real(kind=real64) :: fjul1, fjul2, t_start, t_stop,day1,day2
      integer :: iyr,imn,ida,ihr,mn,iday,iyr1,mn1,iyr2,mn2,unitnum,wout
      real(kind=real32) :: sec,yr1,yr2
      if ((net .eq. 'otr') .or. (net .eq. 'doy')) then
        print*,' start and stop times; year day_of_yr year day_of_yr'
        print*,'  "day_of_yr" may be decimal day'
        read(unitnum,*)yr1,fjul1,yr2,fjul2
        if (wout .eq. 1 ) write(25,*)yr1,fjul1,yr2,fjul2
        call doy2num(int(yr1),int(fjul1),iday)
! days since 1800
        t_start=dble(iday)+fjul1-int(fjul1)

        call doy2num(int(yr2),int(fjul2),iday)
! days since 1800
        t_stop=dble(iday)+fjul2-int(fjul2)
      end if 
      if (net .eq. 'otx' ) then
        print*,' start and stop times; year mn da  year mn da'
        print*,'   "da" may be a decimal day'
        read(unitnum,*) yr1,mn1,day1,yr2,mn2,day2
        if (wout .eq. 1 ) write(25,*) yr1,mn1,day1,yr2,mn2,day2
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        call date2numJL(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'otd' ) then
        print*,' start and stop times; YearMmnDa  YearMmnDay'
        print*,'   "da" may be a decimal day'
        read(unitnum,*) day1, day2
        if (wout .eq. 1 ) write(25,*) day1, day2
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        yr2=float(int(day2/10000.))
        mn2=int((day2-10000*yr2)/100.)
        day2=day2-10000.*yr2- 100.0*mn2
        call date2numJL(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'mjd' ) then
        print*,' start and stop times; MJD  MJD'
        print*,'   May be decimal days'
        read(unitnum,*) day1, day2
        if (wout .eq. 1 ) write(25,*) day1, day2
        t_start=day1-36933.0
        t_stop=day2-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*," start and stop times in GMT time format"
        print*,"   two entries required"
        read(unitnum,*) string1, string2
        if (wout .eq. 1 ) write(25,*) string1, string2
        read(string1,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec
        call date2numJL(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)

        read(string2,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec
        call date2numJL(iyr,imn,ida,iday)
        t_stop=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)
      end if
      return
      end
      subroutine GetTime2(net,n,t_start,t_stop,unitnum,wout)
!
!   read time interval to estimate a RATE change
!
!  net is input for the time format
!  n is an index of the rate change
!  t_start and _stop are output: time span, in days since 1960, to analyze data
!  wout out time on unit 25
!
      use iso_fortran_env
      character(len=3) :: net
      character(len=80) :: string1,string2
      real(kind=real64) :: fjul1, fjul2, t_start, t_stop,day1,day2
      integer :: n, iyr,imn,ida,ihr,mn,iday,iyr1,mn1,iyr2,mn2,unitnum,wout
      real(kind=real32) :: sec,yr1,yr2
      if ((net .eq. 'otr') .or. (net .eq. 'doy')) then
        print*,'rate change interal',n,' year day_of_yr year day_of_yr'
        print*,'  "day_of_yr" may be decimal day'
        read(unitnum,*)yr1,fjul1,yr2,fjul2
        if (wout .eq. 1 ) write(25,*)yr1,fjul1,yr2,fjul2
        call doy2num(int(yr1),int(fjul1),iday)
! days since 1800
        t_start=dble(iday)+fjul1-int(fjul1)

        call doy2num(int(yr2),int(fjul2),iday)
! days since 1960
        t_stop=dble(iday)+fjul2-int(fjul2)
      end if 
      if (net .eq. 'otx' ) then
        print*,'  rate change interal',n,' year mn da  year mn da'
        print*,'   "da" may be a decimal day'
        read(unitnum,*) yr1,mn1,day1,yr2,mn2,day2
        if (wout .eq. 1 ) write(25,*) yr1,mn1,day1,yr2,mn2,day2
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        call date2numJL(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'otd' ) then
        print*,' rate change interal',n,' YearMmnDa  YearMmnDay'
        print*,'   "da" may be a decimal day'
        read(unitnum,*) day1, day2
        if (wout .eq. 1 ) write(25,*) day1, day2
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
        yr2=float(int(day2/10000.))
        mn2=int((day2-10000*yr2)/100.)
        day2=day2-10000.*yr2- 100.0*mn2
        call date2numJL(int(yr2),int(mn2),int(day2),iday)
        t_stop=dble(iday)+day2-int(day2)
      end if
      if (net .eq. 'mjd' ) then
        print*,'  rate change interal',n,' MJD  MJD'
        print*,'   May be decimal days'
        read(unitnum,*) day1, day2
        if (wout .eq. 1 ) write(25,*) day1, day2
        t_start=day1-36933.0
        t_stop=day2-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*,'  rate change interal',n,' in GMT time format'
        print*,'   two entries required'
        read(unitnum,*) string1, string2
        if (wout .eq. 1 ) write(25,*) string1, string2
        read(string1,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec

        call date2numJL(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)

        read(string2,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec
        call date2numJL(iyr,imn,ida,iday)
        t_stop=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)
      end if
      return
      end
      subroutine GetTime3(net,n,t_start,unitnum,wout)
!
!   read time to estimate an OFFSETs
!
!  net is input for the time format
!  n is an index of the rate change
!  t_start day since 1800
      use iso_fortran_env
      character(len=3) :: net
      character(len=80) :: string1
      real(kind=real64) :: fjul1, t_start,day1
      integer :: n, iyr,imn,ida,ihr,mn,iday,iyr1,mn1,iyr2,mn2,unitnum,wout
      real(kind=real32) :: sec,yr1,yr2
      if ((net .eq. 'otr') .or. (net .eq. 'doy')) then
        print*,' Offset time',n,' year day_of_yr '
        print*,'  "day_of_yr" may be decimal day'
        read(unitnum,*)yr1,fjul1
        if (wout .eq. 1 ) write(25,*)yr1,fjul1
        call doy2num(int(yr1),int(fjul1),iday)
! days since 1800
        t_start=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        print*,'  Offset time',n,' year mn da  '
        print*,'   "da" may be a decimal day'
        read(unitnum,*) yr1,mn1,day1
        if (wout .eq. 1 ) write(25,*) yr1,mn1,day1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        print*,' Offset Time',n,' YearMmnDa  '
        print*,'   "da" may be a decimal day'
        read(unitnum,*) day1
        if (wout .eq. 1 ) write(25,*) day1
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        print*,'  Offset time',n,' MJD '
        print*,'   May be decimal days'
        read(unitnum,*) day1
        if (wout .eq. 1 ) write(25,*) day1
        t_start=day1-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*,'  Offset time',n,' in GMT time format'

        read(unitnum,*) string1
        if (wout .eq. 1 ) write(25,*) string1
        read(string1,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec
        call date2numJL(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)
      end if
      return
      end
      subroutine GetTime4(net,n,t_start,unitnum,wout)
!c
!c   read time  to estimate an exponential/Omori trend

!  net is input for the time format
!  n is an index of the rate change
!  t_start day since 1800
!
      use iso_fortran_env
      character(len=3) :: net
      character(len=80) :: string1
      real(kind=real64) :: fjul1, t_start,day1
      integer :: n, iyr,imn,ida,ihr,mn,iday,iyr1,mn1,iyr2,mn2,unitnum,wout
      real(kind=real32) :: sec,yr1,yr2
      if ((net .eq. 'otr') .or. (net .eq. 'doy')) then
        print*,' Exponential time',n,' year day_of_yr '
        print*,'  "day_of_yr" may be decimal day'
        read(unitnum,*)yr1,fjul1
        if (wout .eq. 1 ) write(25,*)yr1,fjul1
        call doy2num(int(yr1),int(fjul1),iday)
! days since 1800
        t_start=dble(iday)+fjul1-int(fjul1)

      end if 
      if (net .eq. 'otx' ) then
        print*,' Exponential time',n,' year mn da  '
        print*,'   "da" may be a decimal day'
        read(unitnum,*) yr1,mn1,day1
        if (wout .eq. 1 ) write(25,*) yr1,mn1,day1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'otd' ) then
        print*,' Exponential Time',n,' YearMmnDa  '
        print*,'   "da" may be a decimal day'
        read(unitnum,*) day1
        if (wout .eq. 1 ) write(25,*) day1
        yr1=float(int(day1/10000.))
        mn1=int((day1-10000*yr1)/100.)
        day1=day1-10000.*yr1- 100.0*mn1
        call date2numJL(int(yr1),int(mn1),int(day1),iday)
        t_start=dble(iday)+day1-int(day1)
      end if
      if (net .eq. 'mjd' ) then
        print*,' Exponential time',n,' MJD '
        print*,'   May be decimal days'
        read(unitnum,*) day1
        if (wout .eq. 1 ) write(25,*) day1
        t_start=day1-36933.0
      end if
      if (net .eq. 'gmt' ) then
        print*,' Exponential time',n,' in GMT time format'

        read(unitnum,*) string1
        if (wout .eq. 1 ) write(25,*) string1
        read(string1,fmt="(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)")iyr,imn,ida,ihr,mn,sec
        call date2numJL(iyr,imn,ida,iday)
        t_start=dble(iday)+dble(ihr)/24.0 + dble(mn)/(24.0*60.0) + dble(sec)/(24.0*3600.0)
      end if
      return
      end
      
      subroutine DataOut(LogNum,netd,time,dat1,dat2)
! output time series of data on assigned logical unit, with specified format
!  LogNum  logical number for output (write(LogNum,fmt=))
!  net  --- time format, doy,otr,gmt,mjd,....
!  itime --  time in days since 1800
!  ftime ---  fractional day
!  dat1 and dat2 are the two 'observations' to output
      use iso_fortran_env
      real(kind=real64) :: time,dec_timed
      integer :: itime,LogNum,nyr,mn,idate,ihr1,imn1,isec1,jul
      real(kind=real32) :: dat1,dat2,sec1
      character(len=3) :: netd
        itime=int(time)
        dec_timed=time-int(time)
          if ((netd .eq. 'otr' ) .or. (netd .eq. 'doy' )) then
            call num2doy(itime,nyr,jul)
            write(LogNum,fmt="(1x,i5,1x, f15.9, 1x, f16.3,1x,f8.2)")nyr,float(jul)+dec_timed,dat1,dat2
          end if
          if (netd .eq. 'otx' ) then
            call num2date(itime,nyr,mn,idate)
            write(LogNum,fmt="(1x,2i5,f15.9, 1x,f16.3,1x,f8.2 )")nyr,mn,idate+dec_timed,dat1,dat2
          end if
          if (netd .eq. 'otd' ) then
            call num2date(itime,nyr,mn,idate)
            write(LogNum,fmt="(1x,i4,i2.2,i2.2,f10.9, 1x,f16.3,1x,f8.2)")nyr,mn,idate,dec_timed,dat1,dat2
          end if 
          if (netd .eq. 'gmt' ) then
            call num2date(itime,nyr,mn,idate)
            ihr1=int(24.0*dec_timed)
            imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
            sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
            sec1=sec1*3600.0*24.0
            isec1=int(sec1)
            write(LogNum,fmt="(1x,i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,f2.1, 1x,f16.3,1x,f8.2)") &
             nyr,mn,idate,ihr1,imn1,isec1,sec1,dat1,dat2
          end if
          if (netd .eq. 'mjd' ) then
            write(LogNum,fmt="(1x,f18.9, 1x,f16.3,1x,f8.2)")dble(itime)+dec_timed+36933.d+0,dat1,dat2
          end if 
      end subroutine DataOut
  
end module GetStuff_mod
  
  