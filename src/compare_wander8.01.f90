Program compare_wander8

! compares the RMS wander or drift of real data with that from from syntheti! data
!   which is generated as either
!  uses two power law functions, white noise, and band-pass filtered white noise

!  USED IN COMPANION with est_noise8.x


!  1) program 1st computes the 'wander' for the real data more or less defined
!      by Agnew, (1992) as function of averaging interval
!     
!  2)  then loops through as series of synthetic data having the
!       same sampling scheme and the same noise model as the real data
!       computes the wander  of the synthetic data
!  3) outputs the wander for the real data, the statistic from syntheti! data.
!       The statistics are the confidence interval (68, 90, and 95) based upon ranking the
!        results of simulation. In addition, the chi^2 confidence is estimated

!  Modified April 2004 for differences between est_noise5fd and est_noise5 in computing
!    covariances.

! Modify Nov 2013 for double precision time

!  This is a re-write of compare_wander6  incorporating two different methods to compute
!  the temporal correlations; either by adding several independent noise sources, as
!   in the original est_noise6/compare_wander6, ie quadrature, or adding the filter functions
!    representing different noise and convolving with a single random noise generator,
!    known as 'additive noise' compatable with est_noise7.2x with 'a'.
! November 2024 -- modified to modern Fortran standards

 use iso_fortran_env
 use GetStuff_mod
 use time_mod
 use compare_wander_mod
 use filterfunc_mod, only: frac_diff, band_pass_filt
 implicit none
 
 real(kind=real64) :: timex,t0,dt_sam
 real(kind=real32), allocatable :: d(:),tim_wan(:),H(:),wan_dat(:),ran1(:),ran(:),dat(:),wan_sim(:), save_wan(:,:),temp(:)
 real(kind=real32) :: dist,err,fmax,tmin,tmax,dtt,dtt1,t_test,tlast,tlog,ftemp,ts,var,sd
 character(len=20) :: filename,name
 character(len=132) :: string
 character(len=3) :: net
 character(len=1) :: ModType
 logical :: filexist 
 integer :: iseed, ios,ic,i,idec,iccc,nmod,ix,nt,ntt,nx,nloop,npts,il,j,n50,n68,n90,n95
 open(25,file='wand_in.jrn')
 open(23,file="seed.dat")
 read(23,*,iostat=ios)iseed
   if ( ios .ne. 0) then
     print*, 'program will use its default value'
     print*, ' best to create a seed file before execution'
     print*, 'echo positive_integer > seed.dat'
     iseed=728745902
     rewind(23)
     write(23,*)iseed

  end if
  
  fmiss=99999.   ! missing data symbol
  print*,'Input the data type for  processing'

  print*,' doy or otr=data'
  print*,'       format of year, day of year, res, calc, obs'
  print*,' otd=data'
  print*,'       format of YearMnDa,  res, calc, obs'
  print*,' otx=data'
  print*,'       format of year mo da, res, calc, obs'
  print*,' mjd=data'
  print*,'       Modified Julian day, res, calc, obs'
  print*,'For all of the above, the da is double precision'
  print*,'    which allows decimal days to accomodate seconds'
  print*,' '
  print*,' gmt=data with GMT format for time'
  print*,'        year-mo-daThr:mn:secs.x , res, calc, obs'
  read(5,*)net

 if ((net .ne. 'mjd') .and.  (net .ne. 'doy') .and. (net .ne. 'otr') .and. & 
      (net .ne. 'otd') .and. (net .ne. 'gmt') .and. (net .ne. 'otx')) then

    print*,' code for network not in list'
     print*,'Must be either doy otr  otx otd mjd or gmt'
    stop
  end if
  write(25,fmt="(a3)")net
  print*,' input the file name of data residuals'
  print*,'  format of data yr doy data; 1999. 80.345 8.1'
  read(5,fmt="(a20)") filename
  write(25,fmt="(a20)")filename


  inquire(file=filename,exist=filexist)
  if ( filexist  ) then
    print*,' data file is present'
  else
    print*,' data file ', filename, ' does not exist'
    print*,'  Balling'
    stop
  end if   !! fileexist?
  open (20,file=filename,action='read')
 !!  get number of data
    ic=0
    do 
      read(20,fmt="(a132)",end=10) string
      call GetData(net,string,timex,dist,err)
      ic=ic+1
    end do  !! end short loop to read and count data 
10  continue
  print*,' number of observations',ic 
  allocate(t(ic))
  allocate(d(ic))
! read the data
  rewind(20)
  do i=1,ic
    read(20,fmt="(a132)") string
    call GetData(net,string,timex,dist,err)
    if (ic .eq. 1) t0=timex
    t(i)=timex-t0
    d(i)=dist
  end do
  
  print*,' input name of file with noise parameters'
  read(5,fmt="(a20)")filename
  open(22,file=filename)
  read(22,fmt="(a20,3x,f10.3,3i6,7(3x,f10.3,1x,10x))")name,fmax,idec,iccc,nmod,sig,amp1,exp1,alpha, & 
     amp_bp,amp2,exp2
  close(22) 
  print*,' Noise parameters read from file'
  print*,sig,amp1,exp1,alpha,amp_bp,amp2,exp2  
  write(25,fmt="(a20)")filename
!  get smallest sampling interval
  dt_sam=999999.
  do  i=2,ic
    if (t(i)-t(i-1) .lt. dt_sam) dt_sam=t(i)-t(i-1)
  end do
  print*,' Sampling interval is ',dt_sam,' days'
  

!  change data spacing 
  ix=int((t(ic)-t(1))/dt_sam) + 1
  allocate(dat(ix))
  dat=fmiss
  do i=1,ic
     ix=int((t(i)-t(1))/dt_sam) + 1
     dat(ix)=d(i)
  end do
  print*,' Number of points spanned by data',ix

!  determine the periods; use constant log(time) to
!   increment the periods
  tmin=dt_sam
  tmax=(t(ic)-t(1))
  dtt=0.05
  dtt1=1./dtt
  nt=int(dtt1*(log10(tmax)-log10(2.0*tmin) )) + 1
  print*, 'nt =',nt
  allocate(tim_wan(nt))
  allocate(wan_dat(nt))
  allocate(wan_sim(nt))

  do i=1,nt
    tlog=float(i-1)*dtt+log10(2.0*tmin)
    if (i .eq. 1) then
        tim_wan(i)=10**tlog

        tlast=tim_wan(i)
        ntt=1
    else
        t_test=10**tlog
        if (t_test-tlast .ge. tmin) then
            ntt=ntt+1
            nx=int(t_test/tmin)
            t_test=float(nx)*tmin
            tim_wan(ntt)=t_test
            tlast=t_test
!            print*,ntt,i,tlog,tim_wan(ntt)
          end if
        end if
  end do
  print*,' Number of intervals to compute wander is:', ntt 
  
  call wander(tim_wan,wan_dat,ntt,dat,ix,dt_sam)  

!  revise list of tim_wan for missing data

  nx=ntt
  do  i=1,ntt
    if (wan_dat(i) .eq. fmiss) then
      nx=nx-1
      do  j=i,nx
          tim_wan(j)=tim_wan(j+1)
          wan_dat(j)=wan_dat(j+1)
      end do
    end if
  end do
  if (nx .ne. ntt) then
    print*,' found ',ntt-nx,' periods without estimating wander'
    ntt=nx
  end if
  print*,' Wander**2 of data is'
  do  i=1,ntt
    print*,tim_wan(i),wan_dat(i)
  end do
 
! Simulate data using noise model

!  read in number of simulations

  print*,' Input the number of sets of sythetic data for simulations'
  read(5,*)nloop
  write(25,fmt="(i5,'     # number of simulations ')")nloop
  print*,' Input the noise model construction'
  print*,'  q -- quadrature add -- convolve each noise filter with'
  print*,'        white noise then add'
  print*,'  n -- same as q'
  print*,'  a -- additive -- add filter functions then convolve'
  print*,'         with white noise'
  read(5,fmt="(a1)")ModType
  write(25,fmt="(a1)")Modtype
  
!  determine filter coeficients for 2 power law noise models, and 1 band-pass filter
!     white noise

  print*,' For bandpass filtered noise, input the pass bands;'
  print*,'   flow and fhigh is c/yrs'
  read(5,*)fl,fh
  if (fl .gt. fh) then
    ftemp=fl
    fl=fh
    fh=ftemp
  end if
  print*,'  and input the number of poles'
  read(5,*)np
  write(25,106)fl,fh,np
106   format(2f10.3, '  #  flow fhigh of BP noise c/yr ',/, i5,'   #  number of poles of BP noise')


!  create the filter functions to creat noise
  allocate(wh(ix+100))
  allocate(pl1(ix+100))
  allocate(pl2(ix+100))
  allocate(bp(ix+100))
  allocate(save_wan(nloop,nt))
  allocate(temp(nloop))
  ts=sngl(dt_sam/365.25)
  npts=ix
  wh=0.0d+0
  pl1=0.0d+0
  pl2=0.0d+0
  bp=0.0d+0

  if (sig .ne. 0) wh(1)=1.0d+0   !  white noise
  if (amp1 .ne. 0) call frac_diff(pl1,exp1,alpha,ts,size(pl1))  ! first PL
  if (amp2 .ne. 0) call frac_diff(pl2,exp2,0.0,ts,size(pl2))  ! second PL                                
  if (amp_bp .ne. 0.0) then
      print*,'calling bp',fl,fh
      call band_pass_filt(dt_sam/365.25,fl,fh,np,npts+100,bp,1.0)
      print*,'bp=',(bp(i),i=1,10)
  end if
  do  il=1,nloop
31    continue      !!  redo if random noise simulation is bad
    call MakeNoise(dat,Modtype,ix,dt_sam,iseed) 

!  variance of data
    var=0.0
    do  i=1,npts
      if (dat(i) .ne. fmiss) var=var+dat(i)**2
    end do
    sd=sqrt(var/float(ic))

    print*,'For synthetic data ',il,' Standard deviation is',sd
    iseed=int(dat(1)*100000) 

    call wander(tim_wan,wan_sim,ntt,dat,ix,dt_sam)
                
! save the results
    do i=1,ntt
      save_wan(il,i)=wan_sim(i) 
    end do
  end do  
  
! Done with computing wander of synthetic data

!  for each period, put the wander in sequential order to do statistics
  do j=1,ntt
    do  i=1,nloop
      temp(i)=save_wan(i,j)
    end do
    call chron(temp,nloop)
    do  i=1,nloop
       save_wan(i,j)=temp(i)
    end do
  end do
  n95=int(nloop*0.025 +0.5)
  n90=int(nloop*0.05 + 0.5)
  n68=int(nloop*0.16 + 0.5)
  n50=int(nloop*0.50 + 0.5)
  print*,' '
  print*,'for wander'
  print*,'68% confidence intervals between ',n68,' and ',nloop-n68
  print*,'90% confidence intervals between ',n90,' and ',nloop-n90
  print*,'95% confidence intervals between ',n95,' and ',nloop-n95
  print*,' median is ',n50
  open(88,file='wander.out')


  do j=1,ntt
    write(88,fmt="(f18.9,9f12.4)")tim_wan(j),sqrt(wan_dat(j)), &
      sqrt(save_wan(n68,j)),sqrt(save_wan(nloop-n68,j)), &
      sqrt(save_wan(n90,j)),sqrt(save_wan(nloop-n90,j)), &
      sqrt(save_wan(n95,j)),sqrt(save_wan(nloop-n95,j)), &
      sqrt(save_wan(n50,j))  !,sqrt(wan_var)

  end do
  print*,' '
  print*,' data wander and confidence limits for plotting in "wander.out" '
  print*,' col 1  is time interval in days, col 2 is sqrt of data wander'
  print*,' col 3 and 4 is 68% CI, 5 and 6 is 90%, and 7 and 8 is 95%'
  print*,'  col 9 is the Median, col 10 is the standard deviation'
  print*,'  col 11 is the confidence that the real data deviates from simulated data'
  print*,'  Journal file of input in wand_in.jrn'
  call random_seed()
  call RANDOM_NUMBER(dist)   !! dist is just a declared variable with no intrinsic meaning
  iseed=int(1000000.*dist-500000)
  rewind (23)
  write(23,*)iseed                                           
  deallocate(wh)
  deallocate(pl1)
  deallocate(pl2)
  deallocate(bp)
  deallocate(dat)
  deallocate(tim_wan)
  deallocate(wan_dat)
  deallocate(wan_sim)
  deallocate(save_wan)
  deallocate(temp)
  deallocate(t)
  deallocate(d)
  close(88)
  close(20)
  close(23)
  close(25)
end program compare_wander8