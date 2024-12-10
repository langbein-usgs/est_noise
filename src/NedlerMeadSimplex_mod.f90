module NedlerMeadSimplex_mod


!   Langbein's coding of Nelder and Mead algorithm of function minimization
!    using Simplex.
!   Adapted for use with est_noise codes
!   Note that est_noise want to maximize, but N-M wants to minimize;
!     hence multiply output of funmin by -1.0
!
!   Note that this code may have similarities with AMOEBA from Numerical Recipes,
!    but I have tried to keep close to the original algorithm.

  use iso_fortran_env
  use est_noise_mod
  use mle_mod
  use alloc_mod
  use model_fit_A_mod

  implicit none
  private
  public :: NedlerMead

  contains
  subroutine NedlerMead(p,y,mp,np,NDIM,FTOL,ITER,npole,flow,fhigh,ModType,xfloat,iswt)

    integer :: iter,mp,ndim,np,NMAX=20,ITMAX,i,ihi,ilo,inhi,j,m,n,xfloat(mp),npole,iswt,itest  
    real(kind=real32) :: y(mp),yswap,ybar,y2,ftol,p(mp,np), nave,rtol,sum,swap,alpha,beta,gam, &
      rtolMin,flow,fhigh
    character(len=1) :: ModType
    real(kind=real64) :: ytry,YPRR,ysave,ytry2
    real(kind=real32), allocatable  :: pave(:),pstr(:),pstr2(:)
    allocate(pave(NMAX))
    allocate(pstr(NMAX))
    allocate(pstr2(NMAX))
    alpha=-1.0
    if (ndim .eq. 1) alpha=-1.1
    beta=0.5
    gam=2.0
    ITMAX=75
    iter=0
    rtolMin=999999.0
    nave=float(ndim+1)
1   continue
    
    do  n=1,ndim
      sum=0.
       do  m=1,ndim+1

           sum=sum+p(m,xfloat(n))/nave
        end do
        pave(xfloat(n))=sum
!        print*,' 1st paves',xfloat(n),pave(xfloat(n))
    end do
    
2   continue 

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
      endif
    end do
 
!  measure deviations  and compare with ftol

! do absolute tolerance 
!      rtol=abs(y(ihi)-y(ilo))
!  extremes
  rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi)+y(ilo)))
!      print*,' rtol=',rtol,y(ihi),y(ilo)
!  RMS
    ybar=0.0d+0
    do  i=1,ndim+1
        ybar=ybar+y(i)/nave
    end do
    y2=0.0d+0
    do i=1,ndim+1
        y2=y2+(y(i)-ybar)**2
    end do
    rtol=(sqrt(y2)/nave)/abs(ybar)
!    print*," rtol rms", rtol
    if (rtol .lt. rtolMin) rtolMin=rtol
!     print*,'rtol=',rtol,y(ihi),y(ilo),ihi,ilo
    if (rtol.lt.ftol) then
        yswap=y(1)
        y(1)=y(ilo)
        y(ilo)=yswap
        do  n=1,ndim
          swap=p(1,xfloat(n))
          p(1,xfloat(n))=p(ilo,xfloat(n))
          p(ilo,xfloat(n))=swap
        end do
!        print*,' Converged! iter= ',iter
        return
      end if

       if (iter.ge.ITMAX) then
         PRINT*,  'NedlerMeade exceeding maximum iterations.',ITMAX
         print*,' rtolMin= ',rtolMin,'rtol=',rtol,y(ihi),y(ilo),ihi,ilo
         return
       end if
      iter=iter+1 
!    Do a reflection and compute y*

     call tryIt (ytry,p,y,pave,mp,np,ndim,pstr,xfloat,ihi,alpha,iswt,ModType,npole,flow,fhigh) 
     

     if (ytry .lt. y(ilo)) then  
            call tryIt2 (ytry2,p,y,pave,mp,np,ndim,pstr,pstr2,xfloat,ihi,gam,iswt,ModType,npole,flow,fhigh)

            iter=iter+1

 
            if (ytry2 .lt. y(ilo))  then
!              print*,' real improvement'
               y(ihi)=ytry2
               do  j=1,ndim
                 p(ihi,xfloat(j))=pstr2(xfloat(j))
               end do
            else
!               print*,'  not improvement'
               y(ihi)=ytry
               do  j=1,ndim
                 p(ihi,xfloat(j))=pstr(xfloat(j))
               end do
            end if
            itest=1
     end if   
     if (ytry .gt. y(ilo)) then
!  Box 2
         itest=0
         if ( (ytry .gt. y(ilo)) .and. (ytry .lt. y(ihi)) ) then
!            print*,' ytry between limits'
!   replace yhi with ytry
            y(ihi)=ytry
            do  j=1,ndim
              p(ihi,xfloat(j))=pstr(xfloat(j))
            end do
            itest=1
         else
            itest=0

!   Do a contraction

           call tryIt(ytry,p,y,pave,mp,np,ndim,pstr,xfloat,ihi,beta,iswt,ModType,npole,flow,fhigh)

           iter=iter+1
!  Box 4
           if (ytry .gt. y(ihi) )  then
           do  i=1,ndim+1
             if(i.ne.ilo) then
               do  j=1,ndim
                 pstr(xfloat(j))=0.5*(p(i,xfloat(j))+p(ilo,xfloat(j)))
                 p(i,xfloat(j))=pstr(xfloat(j))
               end do

               call mle(YPRR,pstr(1),pstr(2),pstr(3),pstr(4),pstr(5),npole,flow,fhigh,pstr(6),pstr(7) &
          ,ModType,iswt)

              iter=iter+1

               y(i)=-1.0*YPRR

             end if
           end do
           itest=1
           else
             itest=1
             y(ihi)=ytry
             do  j=1,ndim
               p(ihi,xfloat(j))=pstr(xfloat(j))
             end do
           end if
         end if

       end if
       if (itest .eq. 1) go to 1

    
    deallocate(pave)
    deallocate(pstr)
    deallocate(pstr2)
     
  end subroutine NedlerMead

  subroutine tryIt(out,p,y,pave,mp,np,ndim,ptry,xfloat,ihi,fac,iswt,ModType,npole,flow,fhigh)


      INTEGER :: ihi,mp,ndim,np,NMAX=20,xfloat(mp),j,iswt,npole,i
      REAL(kind=real32) :: fac,p(:,:),pave(:),fac1,fac2,ptry(:),flow,fhigh,y(:)
    
      real(kind=real64) :: YPRR,ytry,out
      character(len=1) :: ModType
      fac1=(1.-fac)


! initializ ptry
      do  i=1,np
        ptry(i)=p(1,i)
      end do

      do  j=1,ndim
        ptry(xfloat(j))=pave(xfloat(j))*fac1+p(ihi,xfloat(j))*fac

!c        print*,j,ptry(j),pave(j),fac1,p(ihi,j),fac
      end do

      call mle(YPRR,ptry(1),ptry(2),ptry(3),ptry(4),ptry(5),npole,flow,fhigh,ptry(6),ptry(7) &
          ,ModType,iswt)


      ytry=-1.0d+0*YPRR

!      tryIt=ytry
        out=ytry
  
      end subroutine tryIt
      
  subroutine tryIt2(out,p,y,pave,mp,np,ndim,ptry,ptry2,xfloat,ihi,fac,iswt,ModType,npole,flow,fhigh)


      INTEGER :: ihi,mp,ndim,np,NMAX=20,xfloat(mp),j,iswt,npole,i
      REAL(kind=real32) :: fac,p(:,:),pave(:),fac1,fac2,ptry(:),ptry2(:),flow,fhigh,y(:)
    
      real(kind=real64) :: YPRR,ytry,out
      character(len=1) :: ModType
      fac1=(1.-fac)

 
     do  i=1,np
        ptry2(i)=p(1,i)
     end do
      do  j=1,ndim

        ptry2(xfloat(j))=pave(xfloat(j))*fac1+ptry(xfloat(j))*fac

      end do

      call mle(YPRR,ptry2(1),ptry2(2),ptry2(3),ptry2(4),ptry2(5),npole,flow,fhigh,ptry2(6),ptry2(7) &
          ,ModType,iswt)




      ytry=-1.0d+0*YPRR

!      tryIt=ytry
        out=ytry
  
      end subroutine tryIt2
      
      




end module NedlerMeadSimplex_mod