module model_fit_A_mod

  use iso_fortran_env
  use alloc_mod
  use est_noise_mod

  implicit none
  private
  
  public :: fmodel_fit_A, model_fit
 
  
  contains
  

  subroutine model_fit(nopt,ModType,chi2)

!   Subroutine check whether exp_fit or omori_fit has any "floats";
!     if so, used downhill simplex to estimate the optimal time constants
!     if not, just runs calcres
!
!  max_data  dimension of data
!  max_parm  dimension of num parameters
!  ic  number of data
!  nmod number of parameters
!  A  design matrix
!  d  data
!  t  time
!  res  residuals
!  covinv  inverse covariance matrix
!  covar  covariance matrix
!  nopt  outputs the model and its error if set to 1
!   x  the model
!   e  the error
!   nexp  number of exponentials both "fix" and "float"
!   bexp  are the time constants for the exponential (floaters are trial time constrants
!   texp  are the T_o of each of the exponentials
!   exp_choice is either "fix" or "float"
!   t_start  time, in days of first pt in time series
!   t_stop    time, in day of last pt in time series
!  chi2  res^t*covinv*res  chi square
!
    integer :: nopt,i,j,k,kk, nsim, iflag,kcnt,nexp,idum,ITER,ibest,ix,izz,kexp,kkexp
    character(len=1) :: ModType 
    real (kind=real32) :: spar(11,10),ran,rand,bexp_orig,ftol,fmin,e1,e2,Amin,Amax,range,xxd
    real (kind=real64), intent(out) :: chi2
    real (kind=real64) :: sumRes, sres(11),tmin,dt,tmax,tuuu
    nexp=n_exp
    
    if (n_exp .ne. 0 ) then
    
        call modify_A
        
    end if 

    call calcres(nopt,ModType,chi2)

    
    sres(1)=chi2
    
!
!  do a second, third, etc solution depending upon the number of "floaters"
!
    if (n_exp .gt. nexp_fix) then
      nsim=0
      kcnt=0
      rewind (23)
      read(23,*)iseed
      call srand(iseed)
      
      tmin=99999.
      do  i=2,ic
        if (i .gt. size(t)) print*,'exceed',i,size(t)
        dt=t(i)-t(i-1)
        if (dt .lt. tmin) tmin=dt
      end do
      tmax=t(ic)-t(1) 
 
      do  k=1,nexp   !  big loop
        if (exp_choice(k) .eq. "float") then
          kcnt=kcnt+1
          iflag=1
          nsim=nsim+1 
          do  kk=1,11
            spar(kk,nsim)=alog10(bexp(k))

          end do
4         continue   !!! repeat if bexp is too big or too small
          ran=rand()
         
          idum=int(abs(1000000*ran))
          bexp_orig=bexp(k)
          if (exp_type(k) .eq. 'm' ) then
            bexp(k)=(dlog10(20.0d+0/365.25d+0)-dlog10(tmin/(5.0d+0*365.25d+0)))*ran &
            + dlog10(tmin/(5.0d+0*365.25d+0))
          end if
          if (exp_type(k) .eq. 'e' ) then
            bexp(k)=(dlog10(tmax*5.d+0/365.25d+0)-dlog10(tmin/(5.d+0*365.25d+0)))*ran &
              + dlog10(tmin/(5.d+0*365.25d+0))
          end if
          bexp(k)= 10.0**bexp(k)

!!  check to see if new exponential is either too big or to small relative to the length of the data
          if ((bexp(k) .lt. tmin/(365.0*10.)) .or. (bexp(k) .gt. 365.*10.*tmax))  then
             bexp(k)=bexp_orig
             go to 4        !!!  bad programming!
          end if
          if (abs(alog10(bexp(k)) -spar(1,nsim)) .lt. 0.05 ) then
            bexp(k)=bexp(k)*0.12*float(kcnt+1)
            write(89,*)' bexp(k) to close to first one; modify'
            write(89,*) ' bexp(1)= ',bexp(1),' bexp(k)=',k,bexp(k)
          end if

       
          call modify_A
         
          call calcres(nopt,ModType,chi2)

          spar(kcnt+1,nsim)=alog10(bexp(k))
          sres(kcnt+1)=chi2
          
          bexp(k)=bexp_orig
          nmod=nmod_orig 
        end if
      end do   ! end big loop on nexp
      rewind(23)
      write(23,*)idum
!!
!!   Use Nedler Meade downhill simplex to compute optimal time constants
!!
      do  k=1,nsim+1
         write(89,*) ' input to Amoeba, k=',k,sres(k),(spar(k,j),10.0**spar(k,j),j=1,nsim)
      end do
      ftol=1.0e-02
!  check to see if we need to iterate to find best fit time constants
      if (iflag .eq. 0) then
        write(13,fmt="(90(3x,e17.10,1x,e17.10))")(x(i),ex(i),i=1,nmod)
        return
      end if
      call NedlerMeadEXP(spar,sres,11,10,nsim,ftol,ITER,ModType)
      do k=1,nsim+1

         write(89,*) ' output from Nedler Meade, k=',k,sres(k), (j,spar(k,j),10.0**spar(k,j),j=1,nsim)
      end do
      write(89,*)" "      
      fmin=1.e+30
      kcnt=0
      ibest=0
      do  i=1,nsim+1

         if (sres(i) .lt. fmin) then
         ibest=i
         fmin=sres(i)
         end if
      end do

      write(89,*)' ibest=',ibest
      izz=1
      do  j=1,nexp
        if (exp_choice(j) .eq. "float" ) then
           bexp(j)=10.0**spar(ibest,izz)
           izz=izz+1
        end if
      end do
      write(89,*)" Best fit ",fmin,(j, bexp(j),j=1,nexp)
      write(89,*) " "

!      nmod=nmod_orig

      call modify_A
      call calcres(0,ModType,chi2)
     
!   put the estimate of exponential into list of x(i)
         ix=0
         do  k=1,nexp
           if (exp_choice(k) .eq. "float") then
           ix=ix+1
           x(nmod+ix)=bexp(k)
           end if
         end do    
       nmod=nmod_orig
      if (nopt .eq. 1 ) then 

        nmod=nmod+n_exp-nexp_fix   
        do j=1,nmod
          Anor(j)=1.0
        end do
        nmod=nmod+0
!  next, add columns associated with the time-constants
        kexp=0
        kkexp=0
        do  kexp=1,nexp
          if (exp_choice(kexp) .eq. "float") then
            kkexp=kkexp+1

            xxd=1.1
            if (bexp(kexp) .gt. 1.0) xxd=0.9

            Amin=1.0e+20
            Amax=-1.0e+20
            nmod=nmod+1
           
            do  j=1,ic
              A(j,nmod)=0.0
              if (exp_type(kexp) .eq. "e") then
                if (t(j) .ge. texp(kexp)-t_start) then

!   use the derrivative here  ( need to be careful about units between days and years for time
  
                  tuuu=(t(j)-(texp(kexp)-t_start))/365.25
                  A(j,nmod)=(-1.0*tuuu/( (bexp(kexp))**2 ))* exp(-1.0*tuuu/(bexp(kexp) ) )
                  A(j,nmod)=A(j,nmod)*x(nmod_orig+kkexp)

                end if

             end if  !!  end 'e'
             if (exp_type(kexp) .eq. "o") then
               if (t(j) .gt. texp(kexp)-t_start) then


                 e1=alog10(abs(bexp(kexp)) +(sngl(t(j)-(texp(kexp)-t_start))/365.25))
                 e2=alog10(abs(bexp(kexp))*xxd + (sngl(t(j)-(texp(kexp)-t_start))/365.25))
                 A(j,nmod)=1.*(e1-e2)/(1.0*(bexp(kexp)-xxd*bexp(kexp)))*x(nmod_orig+kkexp)

                end if
             end if  !! type 'o'
             if (exp_type(kexp) .eq. "m" ) then
                if (t(j) .gt. texp(kexp)-t_start) then

                   tuuu=(t(j)-(texp(kexp)-t_start))/365.25
                   A(j,nmod)=x(nmod_orig+kkexp)*(-1.0)*tuuu/(((bexp(kexp))**2)*(1.0+tuuu/bexp(kexp))*alog(10.0))
!                  print*,tuuu,A(j,nmod)

                 end if
             end if  !! type "m"
             if (A(j,nmod) .lt. Amin) Amin=A(j,nmod)
             if (A(j,nmod) .gt. Amax) Amax=A(j,nmod)

            end do  !! index j=1,ic
            range=Amax-Amin
            Anor(nmod)=1.0/(range)

!   renormalize the nmod column
           do  j=1,ic
              A(j,nmod)=A(j,nmod)*Anor(nmod)
            end do
!     if bexp is 'out of range' zero out the A matrix

           if ((bexp(kexp) .le. tmin/(365.0*10.)) .or. (bexp(kexp) .gt. 5.*tmax/365.))  then
               do  j=1,ic
                 A(j,nmod)=0.0
                end do
            end if
          end if  !!float
        end do  !!  index kexp=1,nexp

        xsave=x
        nopt=1

        call calcres(nopt,ModType,chi2)

        nopt=0
        do i=1,nmod
!          print*,i,xsave(i),x(i),ex(i),Anor(i)
          ex(i)=ex(i)*Anor(i)
          x(i)=xsave(i)
        end do
        nmod=nmod_orig
      end if   !! nopt eq 1 ??
    end if !! n_exp .gt. nexp_fix
  end subroutine model_fit
  
  subroutine fmodel_fit_A(nopt,ModType,parm,out)
!  max_data,max_parm,ic,nmod,A,t,d,res,
!     &  covinv,covar,iswitch,x,e,
!     &  nexp,bexp,texp,exp_choice,exp_type,t_start,t_stop,parm)

!  function called by Downhill simplex; the sum of square residuals are
!    output to AMOEBA_mod
! max_data  dimension of data
!  max_parm  dimension of num parameters
!  i!  number of data
!  nmod number of parameters
!  ax  design matrix
!  d  data
!  t  time
!  res  residuals
!  covinv  inverse covariance matrix
!  covar  covariance matrix
!  iswitch  outputs the model and its error if set to 1
!  x  the model
!  e  the error
!  nexp  number of exponentials both "fix" and "float"
!  bexp  are the time constants for the exponential (floaters are trial time constrants
!  texp  are the T_o of each of the exponentials
!  exp_choice is either "fix" or "float"
!  t_start  time, in days of first pt in time series
!  t_stop    time, in day of last pt in time series
!   parm is the trial time constants for expontials/omori

      real(kind=real64) :: chi2, sumRes,out
      real(kind=real32) :: parm(20)
      character(len=1) ModType
      integer :: i,ix,nopt,nexp

      nexp=n_exp

      ix=0

      do  i=1,nexp
      if (exp_choice(i) .eq. "float") then
         ix=ix+1

         bexp(i)=10.0**parm(ix)

      end if
      end do

      nmod=nmod_orig
      call modify_A
      call calcres(nopt,ModType,chi2)
      nmod=nmod_orig
      nmod=nmod_orig
      out=chi2

      write(89,8901)chi2,(texp(ix), parm(ix),10.0**parm(i),i=1,ix)
        
8901  format('chi2=',f20.3,5(' texp', f10.3,' log(tau)',f10.3,' tau',e10.3))

      end subroutine fmodel_fit_A
      
!!##########################  Stuff for estimating time constants for exponential and omori law
!!!!!!    Nedler Meade Simplex adapted to this optimization

      subroutine NedlerMeadEXP(p,y,MP,NP,NDIM,FTOL,ITER, ModType)

!   Langbein's coding of Nelder and Mead algorithm of function minimization
!    using Simplex.
!   Adapted for use with est_noise codes
!   Note that est_noise want to maximize, but N-M wants to minimize;
!     hence multiply output of funmin by -1.0

!   Note that this code may have similarities with AMOEBA from Numerical Recipes,
!    but I have tried to keep close to the original algorithm.

!  Modified NedlerMead.f such that the function is the time-dependent model where
!    the time constant of the omori or exponential is being estimated

      INTEGER :: iter,mp,ndim,np,NMAX=20,ITMAX,i,ihi,ilo,inhi,j,m,n,itest,iopt=0
      real(kind=real64) :: y(MP),yswap, ybar,y2,tryIt2,ytry,YPRR,ysave,ytry2
      
      real(kind=real32) ftol,p(mp,np), nave,alpha,beta,rtolMin,gam
      
      character(len=1) ModType
      real(kind=real32) :: rtol,sum,swap,pave(20),pstr(20),pstr2(20)

!  Note, my alpha is -1*alpha of nedler and mead --- allows single subroutine Pave*(1-val)+P*val to be evaluated
      alpha=-1.0
!  If there is only one unknown, then nedler mead doesn't iterate... it becomes
!   degenerate; it automatically picks the best of two solutions,
!   so, to avoid this, make alpha not equal to -1.0
      if (ndim .eq. 1) alpha=-1.1
      beta=0.5
      gam=2.0
!  experiment
      alpha=-1.0 - 0.1/float(ndim)
      beta=0.5 + 0.05/float(ndim)
      gam=2.0+ 0.2/float(ndim)
      ITMAX=25
      iter=0
      rtolMin=999999.0
      nave=float(ndim+1)

1     continue
      do  n=1,ndim
        sum=0.
!        print*,(p(m,(n)),m=1,ndim+1)
        do  m=1,ndim+1
          sum=sum+p(m,(n))/nave
        end do
        pave((n))=sum
      write(89,*)' 1st paves',(n),pave((n))
    end do

2     continue 
!       if (iter .gt. 75) stop 
      ilo=1
!  figure out the y(ihi) and y(ilo)
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      end if
      do  i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        end if
      end do
!       write(89,*)' ihi ilo ',ihi,y(ihi),ilo,y(ilo)


!  measure deviations  and compare with ftol

!  do absolute tolerance 
!      rtol=abs(y(ihi)-y(ilo))
!  extremes
      rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi)+y(ilo)))

      if (rtol .lt. rtolMin) rtolMin=rtol

      if (rtol.lt.ftol) then
        yswap=y(1)
        y(1)=y(ilo)
        y(ilo)=yswap
        do  n=1,ndim
          swap=p(1,(n))
          p(1,(n))=p(ilo,(n))
          p(ilo,(n))=swap
        end do
        write(89,*)' Converged! iter= ',iter
        return
      end if
!      if (iter.ge.ITMAX) pause 'ITMAX exceeded in NedlerMeadeEXP'
       if (iter.ge.ITMAX) then
         write(89,*)'NedlerMeadEXP exceeding maximum iterations.',ITMAX
         write(89,*)' rtolMin= ',rtolMin,'rtol=',rtol,y(ihi),y(ilo)
         return
       end if
      iter=iter+1
!    Do a reflection and compute y*
        call tryItEXP(p,y,pave,mp,np,ndim,pstr,ModType,ihi,alpha,ytry)

!   Box 1
       if (ytry .lt. y(ilo)) then
!   improvement

            call tryItEXP2(p,y,pave,mp,np,ndim,pstr,pstr2,ModType,ihi,gam,ytry2)
            iter=iter+1

!   Box 5
 
            if (ytry2 .lt. y(ilo))  then

               y(ihi)=ytry2
               do  j=1,ndim
                 p(ihi,(j))=pstr2((j))
               end do
            else

               y(ihi)=ytry
               do  j=1,ndim
                 p(ihi,(j))=pstr((j))
               end do
            end if
            itest=1
       end if
!  still part of box 1
       if (ytry .gt. y(ilo)) then
!  Box 2
         itest=0
         if ( (ytry .gt. y(ilo)) .and. (ytry .lt. y(ihi)) ) then
!            write(89,*)' ytry between limits'
!   place yhi with ytry
            y(ihi)=ytry
            do j=1,ndim
              p(ihi,(j))=pstr((j))
            end do
            itest=1
         else
            itest=0

!   Do a contraction

           call tryItEXP(p,y,pave,mp,np,ndim,pstr,ModType,ihi,beta,ytry)

           iter=iter+1

!  Box 4
           if (ytry .gt. y(ihi) )  then
           do  i=1,ndim+1
            if(i.ne.ilo)then
              do j=1,ndim
                pstr((j))=0.5*(p(i,(j))+p(ilo,(j)))
                p(i,(j))=pstr((j))
              end do


              iter=iter+1
            call fmodel_fit_A(iopt,ModType,pstr,y(i))

            end if
           end do
           itest=1
           else
             itest=1
             y(ihi)=ytry
             do j=1,ndim
               p(ihi,(j))=pstr((j))
             end do

           end if
         end if

       end if
       if (itest .eq. 1) go to 1

      return
      END
      subroutine tryItEXP(p,y,pave,mp,np,ndim,ptry,ModType,ihi,fac,result)



      INTEGER :: ihi,mp,ndim,np,NMAX=20,j,i,iopt=0
      REAL(kind=real32) :: fac,p(mp,np),pave(np),fac1,fac2,ptry(20)

      character(len=1) ModType

      
      real(kind=real64) :: YPRR,ytry,result,y(mp)
      fac1=(1.-fac)

!  initializ ptry
      do  i=1,np
        ptry(i)=p(1,i)
      end do
      do  j=1,ndim
        ptry((j))=pave((j))*fac1+p(ihi,(j))*fac

      end do

      call fmodel_fit_A(iopt,ModType,ptry,ytry)
      result=ytry
!      return
      end subroutine tryItEXP
      subroutine tryItEXP2(p,y,pave,mp,np,ndim,ptry,ptry2,ModType,ihi,fac,res)


      INTEGER :: ihi,mp,ndim,np,NMAX=20,j,iopt=0
      REAL(kind=real32) :: fac,p(mp,np),pave(np),fac1,fac2,ptry(20),ptry2(20)
      integer :: i
      character(len=1) ModType
      
      real(kind=real64) :: YPRR,ytry,res,y(mp)
      fac1=(1.-fac)
!  initializ ptry
      do  i=1,np
        ptry2(i)=p(1,i)
      end do
      do  j=1,ndim

        ptry2((j))=pave((j))*fac1+ptry((j))*fac
      end do


      call fmodel_fit_A(iopt,ModType,ptry2,res)

      END subroutine tryItEXP2
 
   
end module  model_fit_A_mod