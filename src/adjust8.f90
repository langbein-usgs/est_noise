!   applies offsets, rate changes, and exponential to time series data
!   accommodates different data formats
program adjust8
  use iso_fortran_env
  use GetStuff_mod
  use time_mod
  implicit none
  character(len=16) :: ifile, ofile
  character(len=1) :: func_type,sgn
  character(len=3) :: netd, nett
  character(len=132) :: string  
  real(kind=real32) :: sn, amp,tau,per,camp,samp,freq,off,rate,sec1
  real(kind=real64) :: tref,toff,toff2,timed,dec_timed
  logical :: filexist
  real(kind=real64), allocatable :: time(:)
  real(kind=real32), allocatable :: data(:),er(:)
  integer :: n,i,noff,nexp,nrate,nlog,nsea,nrates,idate,ihr1,imn1,isec1,itime,jul,mn,nyr,unitnum
  print*,' Program performs adjustments to time series data'
  print*,'   works with offsets, rate changes and exponentials'
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
  print*,' doy or otr   format of year day_of_year obs err'

  print*,' otd   format of YearMnDa  obs err'

  print*,' otx   format of year mo da obs err'

  print*,' mjd   Modified Julian day, obs, err'
  print*,'For all of the above, the da and doy are double precision'
  print*,'    which allows decimal days to accomodate seconds'
  print*,' '
  print*,' gmt  GMT format for time'
  print*,'        year-mo-daThr:mn:secs.x  obs  err'
  read(5,fmt="(a3)")netd
  print*,' Input name of time series data file to apply offsets and adjustments (e.g.  cmbb.u)'
  read(5,fmt="(a16)")ifile
  inquire(file=ifile,exist=filexist)
  if ( filexist  ) then
     print*,' data file is present'
  else
     print*,' data file ', ifile, ' does not exist'
     print*,'  Balling'
     stop
  end if
  open(1,file=ifile,action='readwrite')
  print*,' Input the format of times used in adjustment file'
  print*,'    See above listing'
  read(5,fmt="(a3)")nett
  print*,' Input name of file with adjustments (e.g. adj_cmbb.u)'
  read(5,fmt="(a16)") ofile
  inquire(file=ofile,exist=filexist)
  if ( filexist  ) then
     print*,' adjustment file is present'
  else
     print*,' adjustment file ', ofile, ' does not exist'
     print*,'  Balling'
     stop
  end if
  print*,'  Input either d/a to indicate whether adjustments'
  print*,'   will (d)etrend data, or (a)dd a trend to data'
  print*,'   a d will multiply given adjustments by -1.0'
  read(5,fmt="(a1)")sgn
  sn=1.0
  if (sgn .eq. 'd') sn=-1.0
  open(2,file=ofile,action='read')  
  unitnum=2
!  find out number of data
  n=0
  do 
     read(1,fmt="(a80)",end=10) 
     n=n+1
  end do
10 continue
  print*,' Number of data read ', n
  rewind(1)
  allocate(time(n))
  allocate(data(n))
  allocate(er(n))

!
!  read the data
!
  do i=1,n
    read(1,fmt="(a80)")string
    call GetData(netd,string,time(i),data(i),er(i))
  end do
  noff=0
  nexp=0
  nrate=0
  nlog=0
  nsea=0
  nrates=0
!
!  read and act on adjustment lopp
!
  do
    read(2,*,end=50) func_type
        if (func_type .eq. 'o') then
!    deal with offset
          call GetTime3(nett,noff+1,toff,unitnum,0)
          read(2,*)  off
!          print*,' Offset at ',toff
          do i=1,n
            if (time(i) .ge. toff) data(i)=data(i)+sn*off
          end do
          noff=noff+1
        end if
        if (func_type .eq. 'r') then
!  deal with rate change
          call GetTime2(nett,nrate+1,toff,toff2,unitnum,0)
!         print*,' rate change between', toff,' and ',toff2
          read(2,*) rate
          do  i=1,n
            if ((time(i) .ge. toff) .and. (time(i) .lt. toff2))  data(i)=data(i) &
              +sn*rate*(time(i)-toff)/365.25
            if (time(i) .ge. toff2) data(i)=data(i) + sn*rate*(toff2-toff)/365.25
          end do

          nrate=nrate+1
        end if
        if (func_type .eq. 'e') then
!  deal with exponentials
          call GetTime3(nett,nexp+1,toff,unitnum,0)
!         print*,' exponential at ',toff
          read(2,*) tau, amp
          do  i=1,n
            if (time(i) .ge. toff) data(i)=data(i)+sn*amp*(1.0 - exp(-(time(i)-toff)/(tau*365.25)) )
          end do
          nexp=nexp+1
        end if
        if (func_type .eq. 'm') then
!  deal with logrithm
          call GetTime3(nett,nlog+1,toff,unitnum,0)
 
          read(2,*) tau, amp
!         print*,' log function at ',toff
          do i=1,n
            if (time(i) .ge. toff) data(i)=data(i)+sn*amp*dlog10(1.0 + (time(i)-toff)/(tau*365.25))
          end do
          nlog=nlog+1
        end if
        if (func_type .eq. 's') then
          call GetTime3(nett,nsea+1,tref,unitnum,0)

          read(2,*) per,camp,samp
          freq=2.0*3.1415926/per
          do  i=1,n
            data(i)=data(i)+sn*camp*cos(freq*(time(i)-tref)) + sn*samp*sin(freq*(time(i)-tref))
          end do
         nsea=nsea+1
        end if 
        if (func_type .eq. 'R' ) then
!  deal with removing the rate
          read(2,*)rate
          do  i=1,n
            data(i)=data(i)+sn*rate*(time(i)-time(1))/365.25
          end do
          nrates=nrates+1
       end if
           
  end do    !!!  finish loop reading and adjusting data
50 continue
!
!  output adjusted data
!
      rewind(1)
      do  i=1,n
        itime=int(time(i))
        timed=(time(i))-int(time(i))
        dec_timed=timed
        if ((netd .eq. 'otr' ) .or. (netd .eq. 'doy' )) then
          call num2doy(itime,nyr,jul)
          write(1,fmt="(1x,i5,1x, f15.9, 1x, f16.3,1x,f8.2)")nyr,float(jul)+dec_timed,data(i),er(i)
        end if
        if (netd .eq. 'otx' ) then
          call num2date(itime,nyr,mn,idate)
          write(1,fmt="(1x,2i5,f15.9, 1x,f16.3,1x,f8.2 )")nyr,mn,idate+dec_timed,data(i),er(i)
        end if
        if (netd .eq. 'otd' ) then
          call num2date(itime,nyr,mn,idate)
          write(1,fmt="(1x,i4,i2.2,i2.2,f10.9, 1x,f16.3,1x,f8.2)")nyr,mn,idate,dec_timed,data(i),er(i)
        end if 
        if (netd .eq. 'gmt' ) then
          call num2date(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_timed)
          imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
          sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          write(1,fmt="(1x,i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,f2.1, 1x,f16.3,1x,f8.2)") &
             nyr,mn,idate,ihr1,imn1,isec1,sec1,data(i),er(i)
        end if
        if (netd .eq. 'mjd' ) then
          write(1,fmt="(1x,f18.9, 1x,f16.3,1x,f8.2)")dble(itime)+dec_timed+36933.d+0,data(i),er(i)
        end if
      end do 
      print*,' number of sinusoids fixed', nsea
      print*,' number of offset fixed', noff
      print*,' number of ratechanges fixed', nrate
      print*,' Number of exponentials fixed', nexp
      print*,' Number of secular rates fixed', nrates
      print*,' Number of logrithm/omori fixed', nlog
      print*,' Adjusted data in file ', ifile
  deallocate(time)
  deallocate(data)
  deallocate(er)
  close(2)
  close(1)
end program adjust8
