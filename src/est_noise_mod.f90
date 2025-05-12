module est_noise_mod

!  contains several subroutines require for est_noise program
  use iso_fortran_env
  use time_mod
  use alloc_mod
  use genNoise_mod! use NedlerMeadSimplexEXP_mod
 

  implicit none
  private
  
  public :: chron, decimate, output, modify_A, outResid, calcres

!   and its subroutines;
  contains


  subroutine chron
!  puts data into chronological order
!  code riped-off from numerical recipes piksrt.f
!
!  nd and md are dimensions of t, d, and a
!  t is time   dtime is same but double precision
!  d is data
!  a is the design matrix
!  ic is number of data
!  nmod number of model parameters

    use iso_fortran_env
    implicit none
    real(kind=real64) :: ax(200),tx
    real(kind=real32) :: dx
    integer :: j,i,k

    do  j=2,ic
         tx=t(j)
         dx=d(j)
         do  k=1,nmod
           ax(k)=A(j,k)
         end do

         do  i=j-1,1,-1
!  if time is a duplicate, add 0.00001 days to value (about 1 second)
         if (t(i) .eq. tx) t(i)=t(i)+0.00001
         if (t(i) .lt. tx) go to 10
           t(i+1)=t(i)
    
           d(i+1)=d(i)
           do  k=1,nmod
             A(i+1,k)=A(i,k)
           end do
         end do
         i=0
10       t(i+1)=tx
        
         d(i+1)=dx
         do  k=1,nmod
           A(i+1,k)=ax(k)
         end do
    end do
    return
    end  subroutine
    


  subroutine calcres(iswitch,ModType,chi2)
  use alloc_mod
  use OMP_LIB
!  max_data  dimension of data
!  max_parm  dimension of num parameters
!  ic  number of data
!  nmod number of parameters
!  ax  design matrix
!  d  data
!  dr  residuals
!  covinv  inverse covariance matrix
!  covar  covariance matrix
!  iswitch  outputs the model and its error if set to 1
!   x  the model
!   ex  the error
!   sum Sum of residuals normalized by covariance  dr^t*covinv*dr
!   ModType;  f, c, n
!
    integer :: iswitch,i,j,k,ier,ilen
    integer :: ie,iflag   ! mdmax,maxnmod,
    character(len=1) :: ModType 
    
    real (kind=real64), intent(out) :: chi2
    real (kind=real64) :: wzd(5000),rat,at,resA,chi2adj,minH,tmp(ic+1),tmp1(ic+1)
    real(kind=real32) :: time0,time1


    if (( ModType .eq. 'f' ) ) then
!
!       execute this by "whitening" both the model and data using Finv
!
!        C^-1 = (F^-1)^t * (F^-1)
!        x=A^-1 d = (A^t C^-1 A)^-1 A^t C^-1 d
!                 = [ A^t (F^-1)^t * (F^-1) A ]^-1 A^t (F^-1)^t (F^-1 d)
!        let F^-1 d = dw
!            and (F^-1) A = Aw

!
!       x= [Aw^t Aw - Aw^t B Aw ]^-1 [Aw^t - Aw^B] dw
!    where B= E (E ^t E)^-1 E^t
!    and E = F^-1 M  where M is a matrix composed of indices of missing data
!    Note the B is computed in funmin and passed onto 
!
!    chi^2 =d^t * C^-1 d
!          =dw^t dw -dw^t *B dw

! figure-out the effective length of filter function... when does it go to zero
    minH=maxval(dsqrt(Finv**2))*1.0e-06
    ilen=ic
    do i=1,ic

       if (dsqrt(Finv(ic+1-i)**2) .ge. minH) then
         ilen=ic+1-i
         exit
       end if
    end do
    if (iswitch .eq. 1) ilen=ic
    ilen=ic      !!!  with parallel openmp, shortening the length of the filter has no measurable affect on cpu time

      Aw=A
!$OMP PARALLEL DO
      do j=1,nmod
!            if ( j .eq. nbt+1 ) call convolvT(mdmax,ic,Finv,A(:,j),Aw(:,j))   !! only used if shortening the length of the filter
          call convolvT(mdmax,ilen,Finv,A(:,j),Aw(:,j))
         do  i=1,ic
            Awt(j,i)=Aw(i,j)
         end do
       end do
!$OMP END PARALLEL DO

!     whiten data
      do  i=1,ic
        din(i)=d(i)
      end do
      tmp=d
      call convolvT(mdmax,ic,Finv,tmp,din)
      do i=1,ic
        dw(i)=din(i)
      end do
 
!    Aw^t * Aw
    
 
        call dgemm('T','N',nmod,nmod,ic,1.0d+0,Aw,mdmax,Aw,mdmax, 0.0d+0,a2,maxnmod)      

        if (nmiss .gt. 0) then
!    Aw^t B Aw
!    a3 = Aw^t B   !!  takes time..

          call dgemm('T','N',nmod,ic,ic,1.0d+0,Aw,mdmax,B,mdmax,0.0d+0,a3,maxnmod)


! a2= (Aw^t B) Aw= a2 - a3 Aw


         call dgemm('N','N',nmod,nmod,ic,-1.0d+0,a3,maxnmod,Aw,mdmax, 1.0d+0,a2,maxnmod)

!    Awt = Aw^t - Aw^t B = Aw^t - a3
!$OMP PARALLEL
!$OMP DO
           do  i=1,nmod
           do  j=1,ic
             Awt(i,j)=Awt(i,j)-a3(i,j)
           end do
           end do
!$OMP END DO
!$OMP END PARALLEL


         end if

! invert Aw^t * Aw


      call dsyev('V','L',nmod,a2,maxnmod,evd,wzd,5000,ier)

      if (iswitch .eq. 1) print*,' Eigenvalues ',(i,evd(i),i=1,nmod)

! examine the eigenvalues
      ie=1
      iflag=0
      do i=nmod-1,1,-1
          if (iflag .eq. 1) go to 212
          rat=evd(i)/evd(nmod)
          rat=abs(rat)
          if (rat .lt. 1.0e-09) go to 213
          ie=ie+1
          go to 212
213       iflag =1
212       continue
     end do
     if (iswitch .eq. 1) print*,' Using ',ie,' out of', nmod,' eigenvalues'
     if (ie .ne. nmod)  print*,' Using ',ie,' out of', nmod,' eigenvalues'

!  construct inverse
     
     if (ier .ne. 0) then

          print*,' eigenvalues ',(evd(i),i=1,nmod)
          print*,' ier= ',ier
          print*,'  if ier=0; success'
          print*,'        <0the i-th argument had an illegal value'
          print*,'        >0the algorithm  failed  to  converge'
          print*,'      do a man ssyev'
      end if
      do i=1,nmod
        do j=1,nmod
          at=0.0
          do  k=nmod,nmod-ie+1,-1
             at=at+a2(i,k)*a2(j,k)/(1.0*evd(k))
          end do
          ai(i,j)=at
        end do
      end do

      

! compute the inverse matrx 

      call dgemm('N','N',nmod,ic,nmod,1.0d+0,ai,maxnmod,Awt,maxnmod,0.0d+0,ainv,maxnmod)
      do  i=1,nmod

        x(i)=0.
        ex(i)=0.
        do  j=1,ic
           x(i)=x(i)+ainv(i,j)*dw(j)
           
        end do

      end do

      dr=matmul(A,x)
      
!     whiten residuals

      tmp=dr
      call convolvT(ic,ic,Finv,tmp,dr)

      chi2=0.0
      do  i=1,ic

          chi2=chi2+(dr(i)-dw(i))**2
       end do
!
!      chi2adj= (dout-dw)^t * B * (dout-dw)
       

        din=matmul(B,dr)-matmul(B,dw)

        chi2adj=0.0
        do  i=1,ic

          chi2adj=chi2adj+ (dr(i)-dw(i))*din(i)

        end do

        chi2=chi2-chi2adj    
       
         
    
     
    else      !!!  ModType  f
!
!      execute the following using standard inversion methods
!
!     compute At*covinv
!
      if (nmod .gt. 0 ) then
      
      call dgemm('T','N',nmod,ic,ic,1.0d+0,A,mdmax,covinv,mdmax,0.0d+0,a1,maxnmod)
!    compute At*covinv*A   A1*A

      call  dgemm('N','N',nmod,nmod,ic,1.0d+0,a1,maxnmod,A,mdmax,0.0d+0,a2,maxnmod)

!    invert At*covinv*A

      call dsyev('V','L',nmod,a2,maxnmod,evd,wzd,5000,ier)
    

!    examine the eigenvalues
      ie=1
      iflag=0
      do  i=nmod-1,1,-1
        if (iflag .ne. 1) then
          rat=evd(i)/evd(nmod)
          rat=abs(rat)
          if (rat .lt. 1.0e-9) then
            iflag=1
          else
            ie=ie+1
          end if
        end if
      end do
      if (iswitch .eq. 1) print*,' Eigenvalues ',(i,evd(i),i=1,nmod)
      if (iswitch .eq. 1) print*,' using ',ie,' out of ',nmod,' eigenvalues'
! examine eigenvalue
      if (ier .ne. 0) then

        print*,' eigenvalues ',(evd(i),i=1,nmod)
        print*,' ier= ',ier
        print*,'  if ier=0; success'
        print*,'        <0the i-th argument had an illegal value'
        print*,'        >0the algorithm  failed  to  converge'
        print*,'      do a man ssyev'
      end if
      do  i=1,nmod
      do  j=1,nmod
        at=0.0d+0
        do  k=nmod,nmod-ie+1,-1
          at=at+a2(i,k)*a2(j,k)/evd(k)
        end do
        ai(i,j)=at
      end do
      end do

      
!    compute the inverse matrx 

      call dgemm('N','N',nmod,ic,nmod,1.0d+0,ai,maxnmod,a1,maxnmod,0.0d+0,ainv,maxnmod)

!  compute the model, x
!!!$OMP PARALLEL 
!!!$OMP DO
      do i=1,nmod

        x(i)=0.d+0
        ex(i)=0.d+0
         do  j=1,ic
           x(i)=x(i)+ainv(i,j)*d(j)
         end do
      end do
!!!$OMP END DO
!!!$OMP END PARALLEL
      
!    compute the residuals
      
!!!$OMP PARALLEL
!!!$OMP  DO
      do i=1,ic
      res(i)=0.
        do  j=1,nmod
          res(i)=res(i)+A(i,j)*x(j)
        end do
      res(i)=d(i)-res(i)
      end do
!!!$OMP END DO
!!!$OMP END PARALLEL
     else
      res=d
     end if   !  execute if nmod > 0


!     calculate RMS residuals


      chi2=0.d+0
!!!$OMP PARALLEL 
!!!$OMP DO
      do  i=1,ic
      resA=0.
        do  j=1,ic
          resA=resA+dble(res(i))*covinv(i,j)*dble(res(j))
         end do!         print*,i,res(i)
      chi2=chi2+resA
      end do
!!!$OMP END DO
!!!$OMP END PARALLEL
!      print*,' chi2=',chi2, sqrt(chi2/float(ic-nmod))
    end if   !!  ModType 
    
! ------------End of Inversion
!
!
! compute the error bars

!    

   if (iswitch .eq. 1) then



      do i=1,nmod
!
!   model cov = [A^t * C^-1 A ]^-1  which is ai from above
!

        ex(i)=sqrt(ai(i,i))

      end do
!
!  Output the model Covariance
!
      open(78,file="covar.out")
      rewind (78)
      write(78,319)
319   format("  Model Covariance",/," Listed approximately in same", &
        "   order as output with exception of floating ",  &
        /,"   omori/exponental amplitude and time constant",/)
      do  i=1,nmod
        write(78,fmt="(i3,' std. err= ',e13.4,' cov= ',82e13.4)")i,sqrt(ai(i,i)),  (ai(i,j),j=1,nmod)
    
      end do

!  Output correlation coefficients
      write(78,318)
318   format(//,"  Model correlation coef.",/," Listed approximately ", "in same", &
       "   order as output with exception of floating ", /,"   omori/exponental amplitude and time constant",/)     
      do i=1,nmod
        write(78,fmt="(i3,' std. err= ',e13.4,' Cross correl ',82f8.4)")i,sqrt(ai(i,i)), &
         (ai(i,j)/(sqrt(ai(i,i))*sqrt(ai(j,j))),j=1,nmod)
    
      end do
      close(78)
    end if   !! iswitch
 
!      if ( nopt .eq. 1) then
!     if ( iswitch .eq. 1) then
!        write(13,130)(x(i),ex(i),i=1,nmod)


!130     format(90(3x,e17.10,1x,e17.10))
!      else
!        write(13,130)(x(i),ex(i),i=1,nmod)               ! ,(x(nmod+i),0.0,i=1,nmod)

!      end if
 

    return
  end subroutine calcres

  subroutine modify_A

!  modifies A matrix to include exponential terms
   integer :: nexp,i, j, kexp, iflagsmall
   real (kind=real64) :: tmax, tmin,dtmin
   real (kind=real32) :: ratmax

!
!  for each exp_choice equal to "float", add an additional column to A
!
    nexp=n_exp
    tmax=t(1)
    tmin=t(1)
    dtmin=99999.
    ratmax=5.0
    nmod=nmod_orig
    
    do  i=2,ic
        if (t(i) .ge. tmax) tmax=t(i)
        if (t(i) .lt. tmin) tmin=t(i)
        if (t(i)-t(i-1) .lt. dtmin) dtmin=t(i)-t(i-1)
    end do

  
    do  kexp=1,nexp
        
      if (exp_choice(kexp) .eq. "float") then
        nmod=nmod+1
        write(89,*)kexp,exp_choice(kexp),texp(kexp),bexp(kexp)
!  test for bexp being too big or too small
        iflagsmall=0
        if (bexp(kexp) .lt.  dtmin/(365.25*ratmax) ) then
          bexp(kexp)= dtmin/(365.25*ratmax) 
          iflagsmall=0
        end if
        if (bexp(kexp) .gt. (tmax-tmin)*ratmax/365.25 ) &
             bexp(kexp)= (tmax-tmin)*ratmax/365.25
        if ((bexp(kexp) .gt. (tmax-tmin)*1.0/365.25) .and. &
              (exp_type(kexp) .eq. "m") ) &
             bexp(kexp)= (tmax-tmin)*1.0/365.25
            
          do j=1,ic
            A(j,nmod)=0.0
            if (exp_type(kexp) .eq. "e") then
            if (t(j) .ge. texp(kexp)-t_start) &
               A(j,nmod)=1.0-exp(-(t(j)-(texp(kexp)-t_start)) &
                /(365.25*bexp(kexp)))
            end if
            if (exp_type(kexp) .eq. "o") then
              if ((t(j) .gt. texp(kexp)-t_start) .and. (iflagsmall .eq. 0)) &
                A(j,nmod)=alog10(abs(bexp(kexp)) + &
                  (sngl(t(j)-(texp(kexp)-t_start))/365.25))
              if ((t(j) .gt. texp(kexp)-t_start) .and. iflagsmall .eq. 1) &
                A(j,nmod)=alog10((sngl(t(j)-(texp(kexp)-t_start))/365.25))
!     &             -alog10(abs(bexp(k)))
             end if
            if (exp_type(kexp) .eq. "m") then
              if ((t(j) .gt. texp(kexp)-t_start) .and. (iflagsmall .eq. 0)) &
                A(j,nmod)=alog10( 1.0 + &
                  (sngl(t(j)-(texp(kexp)-t_start))/365.25)/abs(bexp(kexp)) )
              if ((t(j) .gt. texp(kexp)-t_start) .and. (iflagsmall .eq. 1)) &
                A(j,nmod)=alog10((sngl(t(j)-(texp(kexp)-t_start))/365.25)) &
                  -alog10(abs(bexp(kexp)))
       end if  !! when 
 
   end do  ! index on j

   end if
  end do 

  end subroutine modify_A
      
  
  

  subroutine decimate(idec)
!  decimates the data (and the A matrix) into an irregular sampling interval such
!    that it provides some aspects of data's original sampling interval;
!   e.g.  idec=1 will keep two adjacent obs and toss-out the third;
!    idec=2 keeps 2 adjacent obs and toss-out the next 2
!    idec=3 keeps 2 adjacent obs and toss-out the next 3
! idec=1 uses 67% of data
! idec=2 uses 50% of data
! idec=3 uses 40% of data

     
    integer :: idec,i,j,icdec
    real(kind=real32), allocatable :: dtemp(:)
    real(kind=real64), allocatable :: ttemp(:),ax(:,:)
    allocate(dtemp(ic))
    allocate(ttemp(ic))
    allocate(ax(ic,nmod))
    icdec=1

    if (idec .eq. 1) then
       do  i=4,ic-4,3
         dtemp(icdec)=d(i-3)
         ttemp(icdec)=t(i-3)        
         do j=1,nmod
           ax(icdec,j)=A(i-3,j)
         end do
         icdec=icdec+1
         dtemp(icdec)=d(i-2)
         ttemp(icdec)=t(i-2)
         do  j=1,nmod
           ax(icdec,j)=a(i-2,j)
         end do
         icdec=icdec+1
       end do
    end if
    
    if (idec .eq. 2) then
      do  i=7,ic-7,6
         dtemp(icdec)=d(i-6)
         ttemp(icdec)=t(i-6)
         do j=1,nmod
           ax(icdec,j)=a(i-6,j)
         end do
         icdec=icdec+1
         dtemp(icdec)=d(i-5)
         ttemp(icdec)=t(i-5)
         do  j=1,nmod
           ax(icdec,j)=a(i-5,j)
         end do
         icdec=icdec+1
         dtemp(icdec)=d(i-3)
         ttemp(icdec)=t(i-3)
         do j=1,nmod
           ax(icdec,j)=a(i-3,j)
         end do
         icdec=icdec+1
      end do
   end if
   
   if (idec .eq. 3) then
     do  i=11,ic-11,10
       dtemp(icdec)=d(i-10)
       ttemp(icdec)=t(i-10)
       do  j=1,nmod
           ax(icdec,j)=a(i-10,j)
       end do
         icdec=icdec+1
         dtemp(icdec)=d(i-9)
         ttemp(icdec)=t(i-9)
         do j=1,nmod
           ax(icdec,j)=a(i-9,j)
         end do
         icdec=icdec+1
         dtemp(icdec)=d(i-7)
         ttemp(icdec)=t(i-7)
         do  j=1,nmod
           ax(icdec,j)=a(i-7,j)
         end do
         icdec=icdec+1
         dtemp(icdec)=d(i-4)
         ttemp(icdec)=t(i-4)
         do  j=1,nmod
           ax(icdec,j)=a(i-4,j)
         end do
         icdec=icdec+1
      end do
    end if
         
      if (idec .ne. 0) then
      print*,' Using ',icdec-1,' out of ',ic,' values'
      ic=icdec-1
      do  i=1,ic
        d(i)=dtemp(i)
        t(i)=ttemp(i)
        do  j=1,nmod
          a(i,j)=ax(i,j)
        end do
      end do
      end if
      deallocate(ttemp)
      deallocate(dtemp)
      deallocate(ax)
      return
      end subroutine decimate
      
      subroutine output


!  output results of model

     integer :: nyr,mn,idate,ihr1,imn1,isec1,ihr2,imn2,isec2,itref,jul,jul0,jul2
     real(kind=real64) :: dec_time1, dec_time2,dec_time,dec_tref
     real(kind=real32) :: sec1,sec2, e_mag,fmag,sec,sea
     integer :: i,ikk,k,itime,itime1,itime2,mn2,nc,nfix,idate2,nyr0,nyr1,nyr2

!
!  five different output format; otr/doy, otd, otx,mjd, and gmt
!
! Modify to output model parameters compatible with the adjust program
!

      rewind (17)
      nc=0
      sea=0.

!  output reference time
      itref=int(t_start)
      dec_tref=t_start-float(itref)
      call num2doy(itref,nyr0,jul0)
!      print*,tref,itref,nyr0,jul0
      if ((net .eq. 'otr') .or. (net .eq. 'doy')) then
 
        write(6,1700)nyr0,jul0+dec_tref

      end if

      if (net .eq. 'otx' ) then
         call num2date(int(t_start),nyr,mn,idate)
         write(6,1701) nyr,mn,idate+dec_tref
      end if
      if (net .eq. 'otd' ) then
         call num2date(int(t_start),nyr,mn,idate)
         write(6,1702) nyr,mn,idate+dec_tref
      end if
      if (net .eq. 'gmt') then
         call num2date(int(t_start),nyr,mn,idate)
          ihr1=int(24.0*dec_tref)
          imn1=int(24.0*60.0*(dec_tref-dble(ihr1)/24.0))
          sec1=dec_tref-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,1703)nyr,mn,idate,ihr1,imn1,sec1
      end if
      if (net .eq. 'mjd' ) then
         write(6,1704)t_start+36933.d+0
      end if
1700  format("Reference epoch:", i6,f16.9)
1701  format("Reference epoch:", i6,i3,f16.9)
1702  format("Reference epoch:", 2i5,f16.9)
1703  format("Reference epoch:",i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',f02.1)
1704  format("Reference epoch:",f21.9)
      do 64 i=1,nbt
        write(6,201)i,x(i),ex(i)
201     format(' Nomimal value for baseline ',i3,f10.2,' +/- ',f10.2)
64    continue
!
!  Secular rate
!
      nc=nc+nbt
      if (nrate .eq. 1) then
        write(6,202)x(nc+1)/rate_norm,ex(nc+1)/rate_norm
        write(17,1501)x(nc+1)/rate_norm
      end if
202   format(' Rate in units per year ',f16.4,' +/- ',f16.4)
1501  format("R",/,f16.4)
      nc=nc+nrate
!
!  Rate changes
!
      if (n_rat_chg .gt. 0) then
        do 65 k=1,n_rat_chg
        itime1=int(rat_chg1(k))
        dec_time1=rat_chg1(k)-float(itime1)
        itime2=int(rat_chg2(k))
        dec_time2=rat_chg2(k)-float(itime2)
        call num2doy(itime1,nyr,jul)
        call num2doy(itime2,nyr2,jul2)
        write(17,1503)k,nyr,jul+dec_time1,nyr2,jul2+dec_time2,x(nc+k)/rat_chng_norm(k)
        if ((net .eq. 'otr') .or. (net .eq. 'doy')) then

          write(6,203)k,nyr,jul+dec_time1,nyr2,jul2+dec_time2,x(nc+k)/rat_chng_norm(k),ex(nc+k)/rat_chng_norm(k)
        end if
        if (net .eq. 'otx' ) then
          call num2date(itime1,nyr,mn,idate)
          call num2date(itime2,nyr2,mn2,idate2)
          write(6,204)k,nyr,mn,idate+dec_time1, &
           nyr2,mn2,idate2+dec_time2,x(nc+k)/rat_chng_norm(k),ex(nc+k)/rat_chng_norm(k)

        end if
        if (net .eq. 'otd') then
          call num2date(itime1,nyr,mn,idate)
          call num2date(itime2,nyr2,mn2,idate2)
          write(6,2041)k,nyr,mn,idate,dec_time1, &
            nyr2,mn2,idate2,dec_time2,x(nc+k)/rat_chng_norm(k),ex(nc+k)/rat_chng_norm(k)

        end if
        if (net .eq. 'gmt') then

          call num2date(itime1,nyr,mn,idate)
          call num2date(itime2,nyr2,mn2,idate2)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          ihr2=int(24.0*dec_time2)
          imn2=int(1440.d+0*(dec_time2-dble(ihr2)/24.0d+0))
          sec2=dec_time2-dble(ihr2)/24.0d+0-dble(imn2)/(1440.0d+0)
          sec2=sec2*3600.0*24.0
          isec2=int(sec2)
          print*,sec2,isec2
          sec2=sec2-float(isec2)
          write(6,2042)k,nyr,mn,idate,ihr1,imn1,isec1,sec1, &
            nyr2,mn2,idate2,ihr2,imn2,isec2,sec2, &
            x(nc+k)/rat_chng_norm(k),ex(nc+k)/rat_chng_norm(k)

        end if    
        if (net .eq. 'mjd' ) then
          write(6,2043)k,dble(itime1)+dec_time1+36933.d+0, &
            dble(itime2)+dec_time2+36933.d+0,  &
             x(nc+k)/rat_chng_norm(k),ex(nc+k)/rat_chng_norm(k)

        end if
65      continue
      end if
      nc=nc+n_rat_chg
203   format(' Rate change number ',i3,' between ',i5,f10.3, &
         ' and ',i5,f10.3,' is ',f16.4,' +/- ',f16.4)
1503  format('r   # rate change num',i2,/,i5,f10.3,1x,i5,f10.3,/,f16.4)
204   format(' Rate change number ',i3,' between ',2i5,f6.3, &
       ' and ',2i5,f6.3,' is ',f12.4,' +/- ',f10.4)

2041   format(' Rate change number ',i3,' between ',i4,i2.2,i2.2,f4.3, &
       ' and ',i4,i2.2,i2.2,f4.3, ' is ',f12.4,' +/- ',f10.4)

2042   format(' Rate change number ',i3,' between ',i4,'-',i2.2,'-',i2.2,'T', &
        i2.2,':',i2.2,':',i2.2,f2.1,' and ',i4,'-',i2.2,'-',i2.2,'T', &
        i2.2,':',i2.2,':',i2.2,f2.1,' is ',f12.4,' +/- ',f10.4)

2043   format(' Rate change number ',i3,' between ',f16.9,' and ', f16.9,' is ',f12.4,' +/- ',f10.4)
!
!  Sinusoids
!
      if (nper .gt. 0) then
         do 66 k=1,nper
         fmag=x(nc+2*k-1)**2 + x(nc+2*k-0)**2
         fmag=sqrt(fmag)
         e_mag=(x(nc+2*k-1)**2)*(ex(nc+2*k-1)**2) +(x(nc+2*k-0)**2)*(ex(nc+2*k-0)**2)
         e_mag=sqrt(e_mag)/fmag
         write(6,205) per(k),x(nc+2*k-1),ex(nc+2*k-1), x(nc+2*k-0),ex(nc+2*k-0), fmag,e_mag
         write(17,1705)k,nyr0,jul0+dec_tref,per(k), x(nc+2*k-1),x(nc+2*k-0)
         if (k .eq. 1) sea=fmag
66       continue
      end if
      nc=nc+2*nper
1705  format("s   # sinusoid num", i2,/,i5,f16.9,/,f16.9,2f12.4)
205   format(' Period of ',f8.3,' days,  cos amp= ',f14.4,' +/-',f12.4, &
     & '  sin amp= ',f14.4,' +/-',f12.4,'  magnitude= ',  f14.4,' +/-',f12.4)
!
! Offsets
!
      if (noff .gt. 0) then
        do 67 k=1,noff
        itime=int(off(k))
        dec_time1=off(k)-float(itime)
        call num2doy(itime,nyr,jul)
        write(17,1706)k,nyr,jul+dec_time1,x(nc+k)

        if ((net .eq. 'otr') .or. (net .eq. 'doy')) then

          write(6,206)k,nyr,jul+dec_time1,x(nc+k),ex(nc+k)
        end if
        if (net .eq. 'otx') then
          call num2date(itime,nyr,mn,idate)
          write(6,207)k,nyr,mn,idate+dec_time1,x(nc+k),ex(nc+k)
        end if
        if (net .eq. 'otd') then
          call num2date(itime,nyr,mn,idate)
          write(6,2071)k,nyr,mn,idate,dec_time1,x(nc+k),ex(nc+k)
         end if
        if (net .eq. 'gmt') then
          call num2date(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,2072)k,nyr,mn,idate,ihr1,imn1,isec1,sec1, x(nc+k),ex(nc+k)
        end if
        if (net .eq. 'mjd') then
          write(6,2073)k,dble(itime)+dec_time1+36933.d+0,x(nc+k),ex(nc+k)

        end if
67      continue
      end if
206   format(' Offset number ',i5,' at ',i5,f10.3, ' is ',f12.4,' +/- ',f10.4)
207   format(' Offset number ',i3,' at ',i5,i3,f8.3,' is ',f12.4,' +/- ' ,f10.4)
2071  format(' Offset number ',i3,' at ',i4,i2.2,i2.2,f4.3,' is ',f12.4,' +/- ' ,f10.4)
2072   format(' Offset number ',i3,' at ',i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,f2.1,' is ',f12.4,' +/- ',f10.4)
2073   format(' Offset number ',i3,' at ',f16.9,f12.4,' +/- ',f10.4)
1706  format("o  # offset num. ",i2,/,i5,f16.9/,f12.4)
!
!  exponentials and/or Omori
!
      nc=nc+noff
      ikk=n_exp
!      nc=nc+n_file_press
      nfix=0      
      if (n_exp .ne. 0) then
!      do 680 k=1,n_exp
!680   if (exp_choice(k) .eq. "float") nfloat=nfloat+1
!      print*,"nfloat=",nfloat
      do  k=1,n_exp
         if (exp_type(k) .eq. "m" ) print*," Omori function"
         if (exp_choice(k) .ne. "float") then
         
         nfix=nfix+1
         itime=int(texp(k))
         dec_time=texp(k)-float(itime)
        call num2doy(itime,nyr,jul)
        write(17,1708)exp_type(k),k,nyr,jul+dec_time,bexp(k),x(nc+k)

        if ((net .eq. 'otr') .or. (net .eq. 'doy')) then

          write(6,208)k,nyr,jul+dec_time,x(nc+k),ex(nc+k),bexp(k)
        end if
        if (net .eq. 'otx') then
          call num2date(itime,nyr,mn,idate)
          write(6,209)k,nyr,mn,idate+dec_time,x(nc+k),ex(nc+k),bexp(k)
        end if
        if (net .eq. 'otd') then
          call num2date(itime,nyr,mn,idate)
          write(6,210)k,nyr,mn,idate,dec_time,x(nc+k),ex(nc+k),bexp(k)
        end if
        if (net .eq. 'gmt') then
          call num2date(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,211)k,nyr,mn,idate,ihr1,imn1,isec1,sec1, x(nc+k),ex(nc+k),bexp(k)
        end if
        if (net .eq. 'mjd') then
           write(6,212)k,dble(itime)+dec_time1+36933.d+0, x(nc+k),ex(nc+k),bexp(k)
        end if
1708  format(A2,"   #exp or omori num.",i2,/,i5,f16.9,/,f10.5,f10.3)
        else
         ikk=ikk+1
         nc=nc+n_file_press
         itime=int(texp(k))
         dec_time=texp(k)-float(itime)
         
        if ((net .eq. 'otr') .or. (net .eq. 'doy')) then
          call num2doy(itime,nyr,jul)
          write(6,2081)k,nyr,jul+dec_time,x(nc+k),ex(nc+k),bexp(k),ex(nc+0+ikk)
        end if
        if (net .eq. 'otx') then
          call num2date(itime,nyr,mn,idate)
          write(6,2091)k,nyr,mn,idate+dec_time,x(nc+k),ex(nc+k), bexp(k),ex(nc+0+ikk)
        end if
        if (net .eq. 'otd') then
          call num2date(itime,nyr,mn,idate)
          write(6,2101)k,nyr,mn,idate,dec_time,x(nc+k),ex(nc+k), bexp(k),ex(nc+0+ikk)
        end if
        if (net .eq. 'gmt') then
          call num2date(itime,nyr,mn,idate)
          ihr1=int(24.0*dec_time1)
          imn1=int(24.0*60.0*(dec_time1-dble(ihr1)/24.0))
          sec1=dec_time1-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
          sec1=sec1*3600.0*24.0
          isec1=int(sec1)
          sec1=sec1-float(isec1)
          write(6,2111)k,nyr,mn,idate,ihr1,imn1,isec1,sec1, x(nc+k),ex(nc+k),bexp(k),ex(nc+0+ikk)
        end if
        if (net .eq. 'mjd') then
           write(6,2121)k,dble(itime)+dec_time1+36933.d+0,x(nc+k),ex(nc+k),bexp(k),ex(nc+0+ikk)
        end if
         nc=nc-n_file_press

        end if
        end do   ! finish loop on k=1,n_exp
208   format(' Exponential number ',i5,' at ',i5,f10.3,' is ',f10.2,' +/- ',f8.2,'  Time constant is: ', e12.5,' yrs')
!     & ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ',f10.5,' yrs')

209   format(' Exponential number ',i3,' at ',i5,i3,f8.3, ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ', e12.5,' yrs')
210   format(' Exponential number ',i3,' at ',i4,i2.2,i2.2,f4.3,' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs')  
211   format(' Exponential number ',i3,' at ',i4,'-',i2.2,'-',i2.2,'T', i2.2,':',i2.2,':',i2.2,f2.1, &
       ' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs')
212   format(' Exponential number ',i3,' at ',f16.9,f12.4,' +/- ',f10.4,'  Time constant is: ', e12.5,' yrs')
2081  format(' Exponential number ',i5,' at ',i5,f10.3, &
       ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ',e12.5,' +/- ',e12.5,' yrs')
2091  format(' Exponential number ',i3,' at ',i5,i3,f8.3, &
       ' is ',f10.2,' +/- ',f8.2,'  Time constant is: ',e12.5,' +/- ',e12.5,' yrs')    
2101  format(' Exponential number ',i3,' at ',i4,i2.2,i2.2,f4.3, &
        ' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs', &
       ' +/- ',e12.5,' yrs') 
2111  format(' Exponential number ',i3,' at ',i4,'-',i2.2,'-',i2.2,'T', &
        i2.2,':',i2.2,':',i2.2,f2.1,' is ',f10.2,' +/- ' ,f8.2,'  Time constant is: ', e12.5,' yrs', &
         ' +/- ',e12.5,' yrs')     
2121  format(' Exponential number ',i3,' at ',f16.9,f12.4,' +/- ',f10.4,'  Time constant is: ', e12.5,' yrs', &
        ' +/- ',e12.5,' yrs') 
      end if

      nc=nc+nfix
      if (n_file_press .ne. 0) then
        do 679 k=1,n_file_press

!        print*,'nc+k',nc+k,'x(nc+k)',x(nc+k),'aux_norm(k)',aux_norm(k)
        write(6,2108)k,x(nc+k)/aux_norm(k),ex(nc+k)/aux_norm(k)
!        print*,k,x(nc+k),ex(nc+k),aux_norm(k)
679     continue
      end if
2108  format(' response of ', i3, ' input function',f13.5,' +/- ',f12.5)


      return
    
  end subroutine output
  
!######################

 subroutine outResid
!
!  Calculates and Outputs residuals in different formats using flag net
!  A -- design matrix (dimensioned ndim by mdim)
!  d --- data vector
!  x -- the model parameters
!  dtime == time, day number
      integer :: i,j,idate,ihr1,imn1,isec1,itime,jul,nyr,mn
      real(kind=real64) :: dec_timed,time
      real(kind=real32) :: calc_m,res_m,sec1




      do  i=1,ic
        if (A(i,1) .ne. 0.0) then
        calc_m=0.
        time=t(i)
        do j=1,nmod
          calc_m=calc_m+A(i,j)*x(j)
        end do
        res_m=d(i)-calc_m
        itime=int(time+t_start)
        dec_timed=(time)+t_start-int(time+t_start)+1.1564e-09
        if ((net .eq. 'otr' ) .or. (net .eq. 'doy')) then
          call num2doy(itime,nyr,jul)

          write(65,fmt="(1x,f6.0,1x, f15.9, 3(1x,f16.3))")float(nyr),float(jul)+dec_timed,res_m,calc_m,d(i)
        end if
        if (net .eq. 'otx' ) then
           call num2date(itime,nyr,mn,idate)
           write(65,fmt="(1x,2i5,f15.9, 3(1x,f16.3))")nyr,mn,idate+dec_timed,res_m,calc_m,d(i)
         end if
         if (net .eq. 'otd') then
           call num2date(itime,nyr,mn,idate)
           write(65,fmt="(1x,i4,i2.2,i2.2,f10.9, 3(1x,f16.3))")nyr,mn,idate,dec_timed,res_m,calc_m,d(i)
         end if
         if (net .eq. 'gmt') then
         call num2date(itime,nyr,mn,idate)
           ihr1=int(24.0*dec_timed)
           imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
           sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
           sec1=sec1*3600.0*24.0
           isec1=int(sec1)
           sec1=sec1-isec1
           write(65,fmt="(1x,i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,f2.1, 3(1x,f16.3))") &
                    nyr,mn,idate,ihr1,imn1,isec1,sec1, res_m,calc_m,d(i)
         end if
        if (net .eq. 'mjd' ) then
           write(65,fmt="(1x,f18.9, 3(1x,f16.3))")dble(itime)+dec_timed+36933.d+0,res_m,calc_m,d(i)
        end if
        end if  !(A(i,1) .ne 0)
      end do



  end subroutine outResid


end module est_noise_mod
  
