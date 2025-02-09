module mle_mod


  use iso_fortran_env
  use alloc_mod
  use est_noise_mod
  use filterfunc_mod
  use model_fit_A_mod

  implicit none
  private
 
  public :: mle
  
  contains 
  
  subroutine mle(fmle,sig1,sig2,sig3,sig4,sig5,npole,flow,fhigh,sig6,sig7,ModType,iswt)


!  computes maximum likelihood for specified inputs

!  fmle Log maximum likelihod
!  sig1  white-noise
!  sig2  random walk noise
!  sig3  exponent of PSD
!  sig4  time-constant of PSD
!  sig5 is amplitude of Band-pass filtered noise
!    npole is the number of poles for BP filter
!    flow and fhigh is the pass-band in c/yr for BP filter
!  sig6  amplitude of second power law noise
!  sig7  exponent of second power law
!  ts_yrs  sampling interval in years
!  t_year  time in days of data (double precision
!  d   data is sequential order
!  txx  time of data obs in years
!  jmax  number of observations
!  max_data and max_parm  dimensions of matrices
!  ax   design matrix
!  nmod  number of parameters in design matrix
!  ipl_flag_1  a flag; if 1, power law covariance matrix has been computed for
!  iswt  switch for yes/no calc of error bars for model
!  irowmiss  --- indices where there are missing data
!  nmiss -  number of missing data;  note that the A and d matrix/vector are zero for each missing term
    character(len=1), intent(in) :: ModType 
    real(kind=real32), intent(inout) :: sig1,sig2,sig3,sig4,sig5,flow,fhigh,sig6,sig7
    integer, intent(in) :: npole, iswt
    real(kind=real64), intent(out) :: fmle
    real(kind=real32) :: z1,z2
    integer :: md,mmax,nexp,iswitch_verbose,jmax,max_time,i,j,k,kkmax,nerror,jmaxTmp
    real(kind=real32) :: tmle0,fs,Po,pow_db,fny,sig1_orig,sig_min,sig4_orig,siginstr,fc,fo, &
      const1,const2,extend,fnaught,fxlow,covmin,covmax,timex,del_sec,sig1_new
    real(kind=real64) :: small=1.0d-15,t_small,ts_yrs,ttlen,dsum,fsmall=1.0d-30,dett,detmiss, chi2
   integer, allocatable :: iroww(:)  
   integer(kind=int32) :: count, count_rate, count_max, count0 !!  arguments for system clock
!  stuff for ModType=c
    integer :: ix,ixx,krow   
    if (.not. allocated(iroww)) allocate(iroww(kmax))
    

    fmle=0.0
    
    md=size(d)  !!+5
    mdmax=size(t_year) !!+5
    irowOffset=0
    nexp=n_exp
    jmax=ic
    max_time=size(t_year)  !!+5
    ts_yrs=dt_sam/365.25d+0

   
    call system_clock(count, count_rate, count_max)
    count0=count
    iswitch_verbose=0


!   if iswitch_verbose=1----lots of output

3890  continue
    if (jmax .gt. md) then
      print*, 'Number of observations, ',jmax,' greater than',' storage ',md
      stop
    end if
    if (nmod .gt. maxnmod) then
       print*, 'Number of model param., ',nmod,' greater than',' storage ',maxnmod
       stop
    end if
    

!  make sure that everything is positive except  but exponent to power law
    sig1=abs(sig1)
    sig2=abs(sig2)
    sig4=abs(sig4)
    sig5=abs(sig5)
    sig6=abs(sig6)
    sig7=abs(sig7)
    if (sig2 .lt. 0.0001) sig2=0.0001
    if (sig6 .lt. 0.0001) sig6=0.0001
!  limit the first index to expmax
    if (sig3 .gt. expmax) sig3=expmax
    
!
!  figure-out the equivalent PSD parameters....the Po term!
!

    fs=1./ts_yrs
    t_small=ts_yrs
    if (sig2 .ne. 0 ) then

      Po=2.0*(6.28**(-1.0*sig3))/(fs**(1-(sig3/2)))*sig2*sig2

      pow_db=10.0*alog(Po)
    end if


!
!  figure-out lowest value of sig1 detectable
!
!  compute smallest sampling interval

    fny=1/(2*t_small)
    sig1_orig=sig1
    fs=1./ts_yrs
    sig_min=sig2*sqrt(2.0)*((2./3.14159)**(sig3/2.))/(fs**(sig3/4))

    if (sig4 .lt. 0) then
        print*,' '
        print*,' Time constant was ',sig4
        sig4=0.0001
        sig4=2.0*3.14159/(60*t_year(jmax))
        print*,' force time constant to be ',sig4
    end if

    sig4_orig=sig4
      

    siginstr=sig1**2
! test for negative amplitude for band-pass filter noise
    if (sig5 .lt. 0) then
         print*,' force BP-filter coef > 0 ',sig5
         sig5=0.01
    end if
    
    if ( iswitch_verbose .eq. 1 ) then
      print*,' '
      
      print*,' Noise model Type ', ModType
      print*,' Initial parameters'
      print*,'  white noise ',sqrt(siginstr),' mm'

      print*,'  power law amp ',sig2,' units/yr^0.25*sig3'
      print*,'  exponent to PDS ', sig3
      print*,'  cut-off freq for 1-GM ', sig4
      print*,' The amplitude of PSD at 1 c/yr in db is: ',pow_db
      print*,'  Nyquist Freq in c/yr', fny

      if (sig4 .ne. 0) then
        fo=sig4/(2*3.14159265)
        print*,' low frequency break point is: ',fo,' c/yr'
        print*,' amplitdue of PSD at lowest frequency (DC) in db is ',pow_db -10*sig3*alog10(fo)
      end if
      if (sig1 .ne. 0) then

        print*,' amplitude of PSD at highest frequency is db is: ',10*alog10(sig1*sig1/fny)
        fc=(((Po)*fny)/(sig1*sig1))**(1./sig3)
        print*,' cross over freq (c/yr) ',fc
      end if
      print*,' BP-filter amplitude: ',sig5
      print*,' BP-filter parameters; Poles= ',npole,' low and high freqs; ',flow,fhigh,' c/yr'
      if (sig6 .ne. 6) then
        print*,' Amplidude of second power law noise ',sig6, ' units/yr^0.25*sig7'
        print*,' Exponent of second power law noise', sig7

      end if

    end if   !! end verbose

    if (kmax .gt. max_time) then
        print*,' kmax = ',kmax,' exceeds max_time of ',max_time
        stop
    end if
    kkmax=kmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((ModType .eq. 'f' ) .or. (ModType .eq. 'c')) then   !!! for additive noise

       if ((ipl_flag_1 .ne. 1) .or. (iswt .eq. 1))  then
         kkmax=kmax
         irowOffset=0

         if ( (sig4 .gt. 0.0) .and. (ModType .eq. 'c')) then
! normally, create 'ghost' data preceding actual time series for GM
!     for GM noise
!   But, for sig4 roughly equal to the length of the time series, it
!    is close to PL noise and those ghost data aren't needed to compute
!    an accurate covariance matrix
           ttlen=(t_year(jmax)-t_year(1))
           fxlow=1/ttlen
           extend=3.0
!  compare sig4 to ttlen
          fnaught=sig4/6.28
          if (fnaught .lt. fxlow/5.0) then
              sig4=0.0
              kkmax=kmax
            else
             extend=3.0
             kkmax=kmax+ifix(extend*1/(sngl(t_small)*sig4))
           end if

           if (kkmax-kmax .gt. 5000) kkmax=kmax+5000
           
!  what happens when sig4 is very close to zero and requests something
!    that exceeds array size..... put limit on kkmax.
!           if (kkmax .gt. max_time) kkmax=max_time
           irowOffset=kkmax-kmax

           if (allocated(f)) then             
             deallocate(f) 
             allocate(f(kkmax))
           end if   
           if (allocated(filtpl1)) then             
             deallocate(filtpl1) 
             allocate(filtpl1(kkmax))
           end if   
    

         end if
         call frac_diff(filtpl1,sig3,sig4,sngl(t_small),kkmax)
         sig4=sig4_orig
       end if
       const1=(sig2**2)*(t_small**(sig3/2.0)) 
       const1=sqrt(const1)

!  create second power law
      if (ipl_flag_2 .ne. 1) then
        call frac_diff(filtpl2,sig7,0.0,sngl(t_small),kmax)
      end if
      const2=(sig6**2)*(t_small**(sig7/2.0)) 
      const2=sqrt(const2) 

! Add all of the filters
!
       do  i=1,kmax
         f(i)=const1*filtpl1(i) + const2*filtpl2(i)+ sig5*filtbp(i)
        
       end do
       
       if (kkmax .gt. kmax) then
          do i=kmax+1,kkmax
            f(i)=const1*filtpl1(i)
          end do
       end if
       f(1)=f(1)+sig1

       do i=1,kmax
           iroww(i)=i
       end do
    end if   !!!  finish creating filter functions for modtype f and c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((ModType .eq. 'q') .or. (ModType .eq. 'n'))  then
!   create covariance matrix for quadrature noise model


      if (ipl_flag_1 .ne. 1)  then

         call pow_law_cov(t_year,jmax,sig1,sig2,sig3,sig4,covar,mdmax,mdmax,t_small)

      end if
       if (iswitch_verbose .eq. 1) then


         call system_clock(count, count_rate, count_max)

       end if
       z1=0.0
       z2=1.0
      if (ipl_flag_2 .ne. 1)  call pow_law_cov(t_year,jmax,z1,z2,sig7,z1,covarpl2,mdmax,mdmax,t_small)

!` construct covariance do to BP-filtered noise

       if (ibp_flag .ne. 1) then
         if (sig5 .ne. 0) call band_pass_cov(t_year,jmax,t_small,flow,fhigh,npole, mmax,md,1.0,covarbp)
       end if

       if ((ipl_flag_1 .eq. 1) .and. (sig2 .ne. 0)) then
         do  i=1,jmax
           do  j=1,i-1
             covar(i,j)=(sig2**2)*covarpl1(i,j)
           end do
           covar(i,i)=(sig2**2)*covarpl1(i,i)+sig1**2
         end do
       end if
       if (sig5 .ne. 0) then
         do  i=1,jmax
           do  j=1,i
             covar(i,j)=(sig5**2)*covarbp(i,j)+covar(i,j)
           end do
         end do
       end if

       if (sig6 .ne. 0) then
         do  i=1,jmax
         do  j=1,i
           covar(i,j)=(sig6**2)*covarpl2(i,j)+covar(i,j)
         end do
         end do
       end if
        do  i=1,jmax
          do  j=1,i
          covinv(i,j)=covar(i,j)
          covar(j,i)=covar(i,j)
          covinv(j,i)=covar(i,j)
          end do

        end do
        ixx=jmax
!        call system_clock(count, count_rate, count_max)
!         print*,' covariances computed', float(count-count0)/float(count_rate)
    end if    !!!!  end ModType q
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    if ( ModType .eq. 'f' ) then
!
!         invert the filter
!         by inverse convolution
          
          finv(1)=1.0/f(1)

         do i=2,kmax
           dsum=0.d+0
           do  j=1,i-1
              dsum=dsum+(finv(j))*(f(i+1-j))
           end do
           finv(i)=-1.0*(dsum/(f(1)))
!!  Needed --- if finv is zero, the construction of covariance matrix takes too long;
           if (abs(finv(i)) .lt. fsmall) finv(i)=fsmall
         end do
!!  test whether finv is the inverse of f
        if (iswitch_verbose .eq. 1) then
          dsum=0.0d+0
          do i=1,kmax
            fmle=0.0d+0
            do  j=1,i
               fmle=fmle+f(j)*finv(i-j+1)
            end do
            if (i .ne. 1) dsum=fmle**2+dsum
            if (i .eq. 1) dsum=(1.0d+0-fmle)**2 + dsum

          end do
          print*,' Inverse test',sqrt(dsum/(kmax))
          end if



! ----------Used Langbein (2017?) modification of Bos et al 2013 -------
     

!  construct eqn 7 of Bos et al
!   C^-1 - C^-1 F ( F^t C^-1 F )^-1 F^t C^-1
!   where F represents the covainace of missing data (the columns in C for which data are missing)
!    C is covariance matrix for which there is no missing data

!  Revise the above;
!   Let F be the matrix composed of finv (inverse filter)
!   Let M be the matrix composed of indices for missing data
!   C^-1 = F^t F
!   Then, let E=F M
!  Re-write eqn 7 of Bos et al interms of new variables;
!   C^-1 - C^-1 M ( M^t C^-1 M )^-1 M^t C^-1
!  plug in C^-1 = F^t F
!  F^t F - F^t F M [ M^t F^t F M ]^-1 M^t F^t F
!  F^t F - F^t E [ E^t E ]^-1 E^t F

!   Let B = E [ E^t E ]^-1 E^t
!   then
!  F^t F - F^t B F

!  When solving for model, x;
!  x = [A^t C^-1 A ] A^t C^-1 d

!  let Aw= F A
!      dw= F d   Whitening of both model and data;

!  x = [ Aw^t Aw - Aw^t B Aw]^-1 [ Aw^t - Aw^t B ] dw

!  So, in this portion, I will calculate B and the determinant;
!   which is n*log(f(1) + trace [ E^t E ]^-1

       covmax=0.0
       covmin=1.0e+20

!
!    the determinant  --- 
       dett=float(jmax)*dlog(f(1))

      if (iswt .ne. 2) then

        if (nmiss .ne. 0 ) then
        

!  Construct the E matrix; E=F M  where F is the inverse filter matrix and M are columns showing missing data indices

          do i=1,nmiss
            
            do j=1,irowmiss(i)-1
              E(j,i)=0.0

            end do
            do  j=1,kmax-irowmiss(i)+1
               k=irowmiss(i)+j-1
               E(k,i)=finv(j)
            end do
          end do   ! index on i
          
!  compute E^t E



          call dgemm('T','N',nmiss,nmiss,kmax,1.0d+0,E,mdmax,E,mdmax,0.0d+0,AAA,nmiss)

!  Invert AAA=E^t E using cholesky routines
          call dpotrf('L',nmiss,AAA,nmiss,nerror)


!  Compute determinant due to missing data --- works
!   in conjunction with section below that does cholesky decomposition
!    of inverse covariance matrix with time series having no gaps...
!    This is replaced by the approximation above that only uses f(1) to get determinant
          detmiss=0.0d+0
          do  i=1,nmiss
            detmiss=detmiss+dlog(AAA(i,i))
          end do
!          print*,' determinant for missing data', detmiss

          dett=dett+detmiss
         

!    continue with inversion
          call dpotri('L',nmiss,AAA,nmiss,nerror)

!    fill in upper trianguangular part of matrix
          do  i=1,nmiss
            do  j=i+1,nmiss
              AAA(i,j)=AAA(j,i)
            end do
          end do
            
!     
!    compute C^-1*F*(F^t*C^-1*F)^-1
!  compute B = E [ E^t E ]^-1 E^t = E AAA E^t
!
!    E * AAA = AAAA

      
      call dsymm('R','L',kmax,nmiss,1.0d+0,AAA,nmiss,E,mdmax,0.0d+0,AAAA,mdmax)


       
       call dgemm('N','T',kmax,kmax,nmiss,1.0d+0,AAAA,mdmax,E,mdmax,0.0d+0,B,mdmax)
                    
        end if   !! nmiss ne 0

         call model_fit(iswt,ModType,chi2)
   
      end if   !! iswt .n.e 2
      
!  compute maximum likelyhood estimate

!  dett is half the true determainant determinant is

    fmle=(-1.0*dett-0.5*chi2) -0.5*float(jmax-nmiss)*dlog(2.0d+0*3.14159265d+0)


    call system_clock(count, count_rate, count_max)
    del_sec=float(count-count0)/float(count_rate)
 



  end if   !! ModType = f

      
  if (ModType .eq. 'c') then
!  Do cholesky to get inverse
!     call cpu_time(timex)
!     call system_clock(count, count_rate, count_max)
!     print*, 'start cholesky time;', float(count-count0)/float(count_rate)   

!   Revise the data and design matrix by tossing out rows of missing data
         
         ix=1
         ixx=1
         do  i=1,nmiss
           krow=irowmiss(i)
           do  k=ix,krow-1
             dxx(ixx)=d(k)
             corrxx(ixx)=t_year(k)

             do  j=1,nmod
               axx(ixx,j)=A(k,j)
             end do         
             ixx=ixx+1
           end do
          ix=krow+1
        end do
        
        do  k=ix,jmax
          dxx(ixx)=d(k)
          corrxx(ixx)=t_year(k)

           do  j=1,nmod
             axx(ixx,j)=A(k,j)
           end do 
           ixx=ixx+1
         end do

        ixx=ixx-1
      
        do  i=1,ixx

           ix=int((corrxx(i)-corrxx(1))/t_small+0.5)+1
           iroww(i)=ix+irowOffset

        end do

        call makeMatrix(iroww,kmax+irowOffset,f,covar,ixx,0,mdmax,mdmax)
        
 
        do  i=1,ixx
          do  j=1,i
          covinv(i,j)=covar(i,j)
          covar(j,i)=covar(i,j)
          covinv(j,i)=covar(i,j)
          end do

        end do
     end if   !!! done with computing covariance for c
!!!!!!!!!!!!!!!!
     if ((ModType .eq. "c") .or. (ModType .eq. "a") .or. (ModType .eq. "n")) then
!  compute inverse covariance for modtype c and q
        call dpotrf('L',ixx,covinv,mdmax,nerror)

!  bail-out if nerror not equal 0
        if (nerror .ne. 0) then
           print*,' covariance is singular, nerror=',nerror  
!   re jigger noise model so that it might become nonsingular...
           sig1_new=sig1*1.26
           if (sig1 .lt. 0.25) sig1_new=0.40
           print*,' replace amplitude of noise,', sig1,' by ',sig1_new
           sig1=sig1_new
           if (sig4 .ne. sig4_orig) sig4=sig4_orig
           
            go to 3890    !!  jump back and redo the covariance     
         end if
         
        dett=0.0
        do i=1,ixx
          dett=dett+dlog(covinv(i,i))
        end do

        call dpotri('L',ixx,covinv,md,nerror)


!  fill in upper trianguangular part of matrix

        do  i=1,ixx
          do j=i+1,ixx
          covinv(i,j)=covinv(j,i)
          end do
        end do

        jmaxTmp=jmax
        jmax=ixx
 
        call model_fit(iswt,ModType,chi2)
        jmax=jmaxTmp
        call system_clock(count, count_rate, count_max)
        del_sec=float(count-count0)/float(count_rate)
   
      end if  !!!  end for ModType=c
      if (allocated(iroww)) deallocate(iroww)

   
     fmle=(-1.0*dett-0.5*chi2) -0.5*float(jmax-nmiss)*dlog(2.0d+0*3.14159265d+0)
    if (sig4 .ne. sig4_orig) sig4=sig4_orig
    if (sig1 .ne. sig1_orig) sig1=sig1_orig
    if (chi2 .gt. 999999999. ) fmle=-999999999.
!      if (sumres .gt. 999999999. ) sumres=999999999.
    if (fmle .lt. -9999999.) fmle=-9999999.
    if (chi2 .lt. 0.) fmle=-9999999.
    write(6,7491)sig1,sig2,sig3,sig4,sig5,sig6,sig7,dett,chi2,fmle,del_sec
    write(12,7490)sig1,sig2,sig3,sig4,sig5,sig6,sig7,dett,chi2,fmle
    ix=0
    do i=1,n_exp
      if (exp_choice(i) .eq. "float") ix=ix+2
    end do
    if (nmod .gt. 0)    write(13,fmt="(90(3x,e17.10,1x,e17.10))")(x(i),ex(i),i=1,nmod+ix)

7490  format(7f13.4,2x,3f13.3)
7491  format(7f13.4,2x,4f13.3)
!  print out the time constants if 'float'
      if (nexp .ne. 0) write(6,7493)(bexp(k),k=1,nexp)
7493  format( 10(' tau= ',e10.2, '  '))
      



  end subroutine mle
  
!end module mle8.01_mod
end module mle_mod



