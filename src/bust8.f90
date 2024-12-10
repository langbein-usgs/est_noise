program bust8
!  Modify bust_5 to Fortran 90
!
! attempts to remove observations that are either too large or too small
!
!  specify input time series
!
!  compute a running median (not running mean)
!    after specifying the period
!  
!  examine the residuals about the median
!
!   identify observations greater or less than a threshold around the median
!   eliminate those observations
!
!  difference between bust_2 and bust_3 is that new verision allows for 2 different
!     date formats
! 
!  Bust_4 uses  more efficient algorithms to take medians and/or sort data
!
!  Modified bust_4 (now _5) to:
!  1) remove sorting routines from Press et al Numerical recipes; substituting
!      routines written by John Burkardt; https://people.sc.fsu.edu/~jburkardt/f77_src/i4lib/i4lib.f
!  2) Provide more variety in data formats
!  3) Get time subroutines from time.f
!
  use iso_fortran_env
  use GetStuff_mod
  use time_mod
  use bust_mod
  implicit none
  character(len=16) :: ifile, ofile
  character(len=3) :: netd
  character(len=132) :: string  
  logical :: filexist
  real(kind=real64) :: tsam_day,time0,tstart,dec_timed,timed
  real(kind=real64), allocatable :: time(:),day(:)
  real(kind=real32), allocatable :: data(:),er(:),dsrt(:),d(:),r(:),err(:)
  real(kind=real32) :: xmed,per, sum,xx1,xx2,thres0,thres,sec1
  real(kind=real32), parameter :: fmiss=9999999.
  integer :: n,i,nobs,nmax,nwind,n_tot_wind,istart,istop,ix,k,n25,n75,n001,n999,n01,n99,n05,n95,n10,n90
  integer :: nbust,itime,nyr,jul,mn,idate,ihr1,imn1,isec1
  
  print*,'  Input the data format type'

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
  print*,' Input name of time series data file to detect and remove outliers'
  read(5,fmt="(a16)")ifile
  inquire(file=ifile,exist=filexist)
  if ( filexist  ) then
     print*,' data file is present'
  else
     print*,' data file ', ifile, ' does not exist'
     print*,'  Balling'
     stop
  end if
  open(1,file=ifile,action='read')
  
  print*,' Input name of output file' 
  read(5,fmt="(a16)")ofile
  open(2,file=ofile)
  open(3,file="reject.out")
  print*,' Input sampling interval in days'
  print*,'   1--hour = ',1./24.d+0
  print*,'   10--minutes =',10./(24.d+0 * 60.0d+0)
  print*,'   1--minute = ',1.0/(24.d+0 * 60.0d+0)
  print*,'   1--second =', 1.0/(24.d+0 * 3600.0d+0)
  read(5,*)tsam_day
  
!  find out number of data
  nobs=0
  do 
     read(1,fmt="(a80)",end=10) 
     nobs=nobs+1
  end do
10 continue
  print*,' Number of data read ', nobs
  rewind(1)

  allocate(time(nobs))
  allocate(data(nobs))
  allocate(er(nobs))
  allocate(dsrt(nobs))
!
!  read the data
!
  do i=1,nobs
    read(1,fmt="(a132)")string
!    print*,i,string
    call GetData(netd,string,time(i),data(i),er(i))
    dsrt(i)=data(i)
  end do
  nmax=int((time(nobs)-time(1))/tsam_day + 0.5) + 1 
  print*,' Number of possible observations',nmax
  print*,' Number of missing observations',nmax-nobs
!  print*,time(1),data(1)
 ! print*,time(nobs),data(n)

!  remove median and a secular rate
  if (nobs .gt. 0 ) then
     call  i4vec_median (nobs, dsrt, xmed )
   else
      xmed=0.0
   end if
   print*,' the median is:', xmed
   allocate(day(nmax))
   allocate(d(nmax))
   allocate(r(nmax))
   allocate(err(nmax))
   do i=1,nmax
     d(i)=fmiss
     r(i)=fmiss
   end do
   
!  redistribute data and time vectors such that each index corresponds
!    to a possible observation; missing observations are flagged with fmiss
   do i=1,nobs
     n=int((time(i)-time(1))/tsam_day + 0.5) + 1
     d(n)=data(i)
     day(n)=time(i)
     err(n)=er(i)
   end do
   
!  remove the median
   do i=1,nmax
     if (d(i) .ne. fmiss) r(i)=d(i)-xmed
   end do
   
   tstart=0

   print*,' Input the window length in days to compute running median'
   read(5,*)per
   nwind=int(per/tsam_day)
   print*,'  Number of points in window', nwind
   n_tot_wind=int(n/nwind)
   print*,'  Total number of windows is', n_tot_wind+1
    

!  compute medians of windows of data with length period 

   do  k=1,n_tot_wind
     istart=k*nwind-nwind+1
     istop=k*nwind
     ix=0
     do i=istart,istop

           if (r(i) .ne. fmiss) then
           ix=ix+1
           dsrt(ix)=r(i)
           end if
      end do ! loop on i
!        xmed=fmedian(ix,dsrt)
!        print*,istart,istop,ix
      if (ix .gt. 0) then
          call  i4vec_median (ix, dsrt, xmed )
      else
          xmed=0.0
      end if
      print*,istart,istop,ix,xmed
!   remove running median
      do  i=istart,istop
        if (r(i) .ne. fmiss) r(i)=r(i)-xmed
      end do ! loop on i
   end do   ! loop on k  
   
!  same as above but for the last window
    istart=istop+1
    istop=n
    ix=0
    do  i=istart,istop
       if (r(i) .ne. fmiss) then
          ix=ix+1
          dsrt(ix)=r(i)
        end if
    end do  ! loop on i

    if (ix .gt. 0) then
      call  i4vec_median (ix, dsrt, xmed )
    else
      xmed=0.0
    end if
    print*,istart,istop,ix,xmed
    do  i=istart,istop
        if (r(i) .ne. fmiss) r(i)=r(i)-xmed
    end do 
    
!   sort the data to get statistics on data

    sum=0.0
    ix=0
    do i=1,nmax
      if (r(i) .ne. fmiss) then
        ix=ix+1
        dsrt(ix)=r(i)
        sum=sum+r(i)**2
!        write(98,*)ix,dsrt(ix)
      end if
    end do
    print*,"  Number of missing observations is ",nmax-ix
    print*, 'Standard deviation of filtered data is ', sqrt(sum/float(ix))
    call i4vec_frac ( ix, dsrt, 1, xx1 )
    call i4vec_frac ( ix, dsrt, ix, xx2 )
    print*," Extreme values ",xx1,xx2, " distance = ",xx2-xx1
    n75=int(0.125*ix)
    n25=ix-n75+1
    call i4vec_frac ( ix, dsrt, n75, xx1)
    call i4vec_frac ( ix, dsrt, n25, xx2)
    print*," 75% interval",xx1, xx2, " distance= ",xx2-xx1
    
    thres0=xx2-xx1
    n90=int(0.05*ix)
    n10=ix-n90+1
    call i4vec_frac ( ix, dsrt, n90, xx1)
    call i4vec_frac ( ix, dsrt, n10, xx2)
    print*," 90% interval",xx1,xx2, " distance= ",xx2-xx1
    
    n95=int(0.025*ix)      
    n05=ix-n95+1
    if (n95 .lt. 1) n95=1
    if (n05 .gt. ix) n05=ix
    call i4vec_frac ( ix, dsrt, n95, xx1)
    call i4vec_frac ( ix, dsrt, n05, xx2)
    print*," 95% interval",xx1,xx2," distance= ",xx2-xx1
    
    n99=int(0.005*ix)
    n01=ix-n99+1
    if (n99 .lt. 1) n99=1
    if (n01 .gt. ix) n01=ix
    call i4vec_frac ( ix, dsrt, n99, xx1)
    call i4vec_frac ( ix, dsrt, n01, xx2)
    print*," 99% interval",xx1,xx2, " distance= ",xx2-xx1
      
    n999=int(0.0005*ix)
    n001=ix-n999+1
    if (n999 .lt. 1) n999=1
    if (n001 .gt. ix) n001=ix
    call i4vec_frac ( ix, dsrt, n999, xx1)
    call i4vec_frac ( ix, dsrt, n001, xx2)
    print*," 99.9% interval",xx1,xx2," distance= ",xx2-xx1
    
    print*, " Input the threshold used to reject data"
    print*,'    use units of IQR'
    read(5,*)thres
    thres=thres*thres0/2.0
    print*,'  The threshold is ',thres
    nbust=0
    
! output cleaned data

   do i=1,nmax
     if (d(i) .ne. fmiss) then

        if (abs(r(i)) .lt. thres) then
!  good data
          call DataOut(2,netd,day(i),d(i),err(i))
 
        else   !  threshold
!  outliers go here
          nbust=nbust+1
          call DataOut(3,netd,day(i),d(i),err(i))
 
        end if   !! threshold
     end if  !! fmiss test
   end do
   
    print*," Number of data rejected is ", nbust


    print*," Clean data in ", ofile
    print*,"  Outlyers in file reject.out"
   deallocate(day)
   deallocate(d)
   deallocate(r)
   deallocate(err)
  deallocate(time)
  deallocate(data)
  deallocate(er)
  deallocate(dsrt)
  
   close(1)
   close(2)
   close(3)
end program bust8
