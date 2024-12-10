module time_mod
!!  Set of subroutines for date conversion.  This will work for a 400
!!  year interval between 1800 and 2200
!!
!  doy2num(year,doy,iday1800)  converts date (year and day of year) to # of days since Dec 31, 1799
!  num2doy(iday1800,year,doy) converts number of days since Dec 31, 1799 to year and doy
!  date2numJL(year,month,day,iday1800)  convert year, month, day to # of days
!  num2date(iday1800,year,month,day) converts number of days since Dec 31, 1799 to year, month and day
!!
!!  Note, MJD day 0 is 1858 11 17
!!
!  ChangeTsam --- rationalize short time intervals (<1-day) to be integer multiples of hour, minute, or second
!  Frac2hrmnsec --- converts fraction of a day to hours,minute,seconds
!!
  implicit none
  private
  public :: doy2num, num2doy, date2numJL, num2date, ChangeTsam, Frac2hrmnsc
  
  contains 
  
  subroutine doy2num(year,doy,iday1800)
!  converts year and doy to number of days since 1800
  implicit none
  integer, intent(in) :: year,doy
  integer, intent(out) :: iday1800
  integer :: leap
  leap=int((year-1801)/4)
  if (year .ge. 1901) leap=leap-1
  if (year .ge. 2101) leap=leap-1
  iday1800=365*(year - 1800)  + leap + doy
  end subroutine

  subroutine num2doy(iday1800,year,doy)
! converts number of days since 1800 to year and doy
  implicit none
  integer, intent(in) :: iday1800
  integer, intent(out) :: year, doy
  integer, parameter :: maxlen=400
  integer :: leap, yrlist(maxlen),i
  leap=0
  do i=1,maxlen
    leap=int((i-2)/4)
    if (i .ge. 1901-1799 ) leap=leap-1
    if (i .ge. 2101-1799 ) leap=leap-1
    yrlist(i)=365*(i-1) + leap + 1
  end do
  i=1
   do while ( (iday1800 >= yrlist(i) ))
    i=i+1
  end do
  i=i-1
  year=i+1799
  doy=iday1800-yrlist(i)+1
  end subroutine
    
  
  subroutine date2numJL(year,month,day,iday1800)
!  convert year, month, day to number of days since Dec 31, 1799
  implicit none
  integer, intent(in) :: year,month,day
  integer, intent(out) :: iday1800
  integer, dimension(12) :: mnth, lmnth
  integer :: i, leap
  character (len=1) :: leapyn
  mnth = [0,31,59,90,120,151,181,212,243,273,304,334]  ! doy-1 for start of each month

!  revise list for leap years
  do i=1,12
    lmnth(i)=mnth(i)
    if (i > 2) lmnth(i)=mnth(i)+ 1
  end do
! get daynumber for Jan 1 of year;
  leap=int((year-1801)/4)


  if (year .ge. 1901) leap=leap-1
  if (year .ge. 2101) leap=leap-1

  iday1800=365*(year - 1800)  + leap + 0
!  Is the year a leap year
   leapyn='n'
   if (float(year)/4.0 - year/4 == 0 ) leapyn='y'
   if (year == 1800) leapyn='n'
   if (year == 1900) leapyn='n'
   if (year == 2100) leapyn='n'
   
   if (leapyn == 'n') then
     iday1800=mnth(month)+day+iday1800
   else
     iday1800=lmnth(month)+day+iday1800    ! for leap years
   end if
  
  end subroutine
  
  subroutine num2date(iday1800,year,month,day)
! convert day number since dec 31, 1799 to year, month, day
  implicit none
  integer, intent(in) :: iday1800
  integer, intent(out) :: year,month,day
  integer :: leap, i,doy
  integer, dimension(13) :: mnth, lmnth
  character (len=1) :: leapyn
  integer, parameter :: maxlen=400
  integer, dimension(maxlen) :: yrlist
  
  mnth = [0,31,59,90,120,151,181,212,243,273,304,334,365]  ! doy-1 for start of each month

!  revise list for leap years
  do i=1,13
    lmnth(i)=mnth(i)
    if (i >= 3 ) lmnth(i)=mnth(i)+1
  end do
  
!  determine the year -- follow step from num2doy subroutine
  leap=0
  do i=1,maxlen
    leap=int((i-2)/4)
    if (i .ge. 1901-1799 ) leap=leap-1
    if (i .ge. 2101-1799 ) leap=leap-1
    yrlist(i)=365*(i-1) + leap + 1
  end do
!  i=1
!   do while ( (iday1800 >= yrlist(i) ))
!    i=i+1
  do i=1,maxlen
    if (yrlist(i) > iday1800) exit
  end do
  i=i-1
  year=i+1799
  day=iday1800-yrlist(i)+ 1
  
!  Is the year a leap year
   leapyn='n'
   if (float(year)/4.0 - year/4 == 0 ) leapyn='y'
   if (year == 1800) leapyn='n'
   if (year == 1900) leapyn='n'
   if (year == 2100) leapyn='n'
   
!  get the month
!  print*,'day=',day
!  i=1
  if (leapyn == 'n') then
!  i=1
!    do while ( (mnth(i) < day  ))
!      i=i+1
    do i=1,13
     if (mnth(i) >= day) exit
    end do
    month=i-1
    day=day-mnth(month)
  else
!    do while ( (lmnth(i) < day ))
!      i=i+1
    do i=1,13
      if (lmnth(i) >= day) exit
    end do
    month=i-1
    day=day-lmnth(month)
  end if 
  
  end subroutine
  
  subroutine ChangeTsam(dt)
!  rationalize the interval dt to be exact multiples of either day, hour, minutes, or seconds
  use iso_fortran_env
  implicit none
  real(kind=real64), intent(inout) :: dt   !  sampling interval in days
  integer :: hr,mn,sec
  print*,dt
  if (dt .lt. 1.0d+0) then
    if (( dt .lt. 1.0d+0) .and. (dt .gt. 4.16d-2)) then
!!   sampling interval is multiples of 1-hour
     hr=int(dt*24.0d+0)
     dt=hr/24.0d+0
    end if
    if ((dt .le. 4.16d-2 ) .and. (dt .gt. 0.000694)) then
!! sampling interval is multiple of 1-min
      mn=int(dt*1440.d+0)
      dt=mn/((24.0d+0 * 60.d+0))
    end if
    if ((dt .le. 0.000694) .and. (dt .gt. 11.55d-6 ) ) then
!! interval is multiple of 1 second
      sec=int(dt*24d+0 * 3600.d+0)
      dt=sec/(24.d+0 * 3600.d+0)
    end if
  end if
  end subroutine
  
  subroutine Frac2hrmnsc(frac,hr,mn,sc)
  use iso_fortran_env
  implicit none
  real(kind=real64), intent(in) :: frac    ! fraction of one day
  integer, intent(out) :: hr,mn,sc
!  real*4, intent(out) :: sc
  real(kind=real64) :: dfrac
  hr=0
  mn=0
  sc=0.0
  if (frac .lt. 1.0d+0) then
    hr=int(frac*24.d+0 + 0.001)
    dfrac=frac - hr/24.d+0
    mn=int(0.001+dfrac*24.0d+0 * 60.d+0)
    dfrac=dfrac-mn/(24.0d+0 * 60.d+0)
    sc=int(0.001+dfrac*(24.0d+0 * 3600.d+0))
  end if
  end subroutine 
  
end module time_mod
  