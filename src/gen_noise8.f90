program gen_noise8
!  program generates a time series of colored (temporally correlated) noise
  use time_mod
  use genNoise_mod
!  use filterfunc_mod
  use iso_fortran_env
!  use iso_fortran_env
  implicit none
  real(kind=real64), allocatable :: t(:)
  real(kind=real64) :: dt_sam, dec_timed   !! dec_timed is fractional part of a day
  real, allocatable :: d(:)
  real :: t_len
  integer :: iseed, ios, npts, i,itime,tnum0,doy,nyr,mn,idate,hr,min,sec
  character(len=1) :: ModType
  character(len=3) :: net 

! input the seed for random number generator
  open(1,file="seed.dat")
  read(1,*,iostat=ios)iseed
  if ( ios .ne. 0) then
    print*, 'no file with value of seed; create a seed.dat file'
  end if

  
  open(25,file='gen_in.jrn')   ! provide a journal file that records the input needed to generate noise

!  specifying the noise, sampling interval, and length of time series

  print*,'  Program generates power law noise; more to follow'
  print*,' Input the method to create power law noise'
  print*,'  n or q --  noise added from independent sources (quadrature)'
  print*,'  a -- filter functions added than convolved with single noise source "additive"'
  read(5,fmt="(a1)")ModType
  if (ModType .eq. "q") ModType="n"
  write(25,fmt="(a1)")Modtype
  print*,' Input the sampling interval in days'
  print*,'   1= 1 day'
  print*,'   0.0416667= 1 hour'
  print*,'   0.0069444444= 10 minutes'
  print*,'   0.00069444444= 1 minute'
  print*,'   1.157407e-05 = 1 second'
  read(5,*)dt_sam
  write(25,fmt="(e16.7,'     # sampling interval in days')")dt_sam
  print*,' Length of time series in years'
  read(5,*)t_len
  write(25,fmt="(f16.7,'     # time series length in years')")t_len
  npts=int(t_len*365.25/dt_sam)
  print*,' Number of point in time-series', npts 

!!  allocate memory for d and t (noisy-data and time)  
  allocate(t(npts))
  allocate(d(npts))
  
!  deal with round-off error in dt_sam
  call ChangeTsam(dt_sam)
  
!!  initialize t and d
  do i=1,npts
     t(i)=(float(i-1))*dt_sam
     d(i)=0.0
  end do

!! continue with input
  print*,' What time format should the noisy data have?'
  print*,'   doy or otr --  year, day_of_year'
  print*,'   otx --  year Mn Da'
  print*,'   otd --  yearMnDay'
  print*,'   gmt --  GMT compatible time'
  print*,'   mjd --  modified julian day'
  read(5,*)net
  write(25,fmt="(a3)")net

!  generate the noise -- subroutine will ask for noise parameters, too
!   call genNoise(ModType,t,d,npts,max_data,iseed,dt_sam,max_time)
   call genNoise(ModType,t,d,iseed,dt_sam)

!  output a new seed for next time program is used
  rewind(1)
  call random_seed()
  call RANDOM_NUMBER(t_len)   !! t_len is just a declared variable with no intrinsic meaning
  iseed=int(1000000.*t_len)
  write(1,*)iseed
  close(1)
  close(25)
  print*,' Journal file of input in gen_in.jrn'

! output the synthetic data


!   use jan 1, 1960 as the first date
  call doy2num(1960,1,tnum0)

  open(2,file='gen.out')
  open(3,file='noise.dat')

  do i=1,npts
    t(i)=t(i)+tnum0
    write(2,*)d(i)
    itime=int(t(i))
    dec_timed=(t(i))-int(t(i))
      if ((net .eq. 'otr' ) .or. (net .eq. 'doy') ) then
        call num2doy(itime,nyr,doy)
        write(3,fmt="(1x,f6.0,1x, f15.9, 1x, f15.3, ' 1.0')")float(nyr),float(doy)+dec_timed,d(i)
      end if
      if (net .eq. 'otx' ) then
        call num2date(itime,nyr,mn,idate)
        write(3,fmt="(1x,2i5,f15.9, 1x,f15.3, '  1.0')")nyr,mn,idate+dec_timed,d(i)
      end if
      if (net .eq. 'otd' ) then
        call num2date(itime,nyr,mn,idate)
        dec_timed=dec_timed+float(idate)+float(100*mn)+float(10000*nyr)
        write(3,fmt="(1x,f19.9,f15.2,' 1.0')") dec_timed,d(i)
      end if
      if (net .eq. 'mjd' ) then
        write(3,fmt="(1x,f18.9, 1x,f15.3, ' 1.0')")dble(itime)+dec_timed-36933.d+0,d(i)
      end if 
      if (net .eq. 'gmt ') then
         call num2date(itime,nyr,mn,idate)
        call Frac2hrmnsc(dec_timed,hr,min,sec)
        write(3,fmt="(1x,i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2, 1x,f15.3, ' 1.0')") &
          nyr,mn,idate,hr,min,sec,d(i)
      end if  
  end do
  close(3)
  close(2)
  close(1)
  deallocate(t)
  deallocate(d)
  print*,' Sequence of correlated data in gen.out; with no time stamps'
  print*,' Simulated noisy data in noise.dat'
  

end program


 
