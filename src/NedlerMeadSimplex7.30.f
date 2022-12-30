      subroutine NedlerMead(p,y,MP,NP,NDIM,FTOL,ITER,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat)
c
c   Langbein's coding of Nelder and Mead algorithm of function minimization
c    using Simplex.
c   Adapted for use with est_noise codes
c   Note that est_noise want to maximize, but N-M wants to minimize;
c     hence multiply output of funmin by -1.0
c
c   Note that this code may have similarities with AMOEBA from Numerical Recipes,
c    but I have tried to keep close to the original algorithm.
c
      dimension d(max_data),res(max_data)
      double precision t_year(max_data),dt_sam,y(mp),yswap,tryIt2,
     &   ybar,y2
      double precision A(max_data,max_mod),tryIt,ytry,YPRR,ysave,ytry2
      INTEGER iter,mp,ndim,np,NMAX,ITMAX,ifloat(mp)
      REAL ftol,p(mp,np), nave
      PARAMETER (NMAX=20)
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,pave(NMAX),pstr(NMAX),pstr2(NMAX)
c  Note, my alpha is -1*alpha of nedler and mead --- allows single subroutine Pave*(1-val)+P*val to be evaluated
      alpha=-1.0
c If there is only one unknown, then nedler mead doesn't iterate... it becomes
c   degenerate; it automatically picks the best of two solutions,
c   so, to avoid this, make alpha not equal to -1.0
      if (ndim .eq. 1) alpha=-1.1
      beta=0.5
      gam=2.0
      ITMAX=75
      iter=0
      rtolMin=999999.0
      nave=float(ndim+1)
c      print*,'simplex', MP,NP,NDIM,FTOL,ITER
1     do 12 n=1,ndim
        sum=0.
c        print*,n,(p(m,ifloat(n)),m=1,ndim+1)
        do 11 m=1,ndim+1
          sum=sum+p(m,ifloat(n))/nave
11      continue
        pave(ifloat(n))=sum
c        print*,' 1st paves',ifloat(n),pave(ifloat(n))
12    continue

2     continue 
c       if (iter .gt. 75) stop 
      ilo=1
c  figure out the y(ihi) and y(ilo)
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
c      print*,' ihi ilo ',ihi,y(ihi),ilo,y(ilo)


c  measure deviations  and compare with ftol

c do absolute tolerance 
c      rtol=abs(y(ihi)-y(ilo))
c  extremes
      rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi)+y(ilo)))
c      print*,' rtol=',rtol,y(ihi),y(ilo)
c  RMS
      ybar=0.0d+0
      do 51 i=1,ndim+1
        ybar=ybar+y(i)/nave
51    continue
      y2=0.0d+0
      do 52 i=1,ndim+1
        y2=y2+(y(i)-ybar)**2
52    continue
      rtol=(sqrt(y2)/nave)/abs(ybar)
c      print*," rtol rms", rtol
      if (rtol .lt. rtolMin) rtolMin=rtol
cs      print*,'rtol=',rtol,y(ihi),y(ilo),ihi,ilo
      if (rtol.lt.ftol) then
        yswap=y(1)
        y(1)=y(ilo)
        y(ilo)=yswap
        do 14 n=1,ndim
          swap=p(1,ifloat(n))
          p(1,ifloat(n))=p(ilo,ifloat(n))
          p(ilo,ifloat(n))=swap
14      continue
c        print*,' Converged! iter= ',iter
        return
      endif
c      if (iter.ge.ITMAX) pause 'ITMAX exceeded in NedlerMeade'
       if (iter.ge.ITMAX) then
         PRINT*,  'NedlerMeade exceeding maximum iterations.',ITMAX
         print*,' rtolMin= ',rtolMin,'rtol=',rtol,y(ihi),y(ilo),ihi,ilo
         return
       end if
      iter=iter+1
c    Do a reflection and compute y*
        ytry=tryIt(p,y,pave,mp,np,ndim,pstr,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat,ihi,alpha)


c   Box 1
       if (ytry .lt. y(ilo)) then
c   improvement
c         print*,' this is an improvement'
c      This is 'expansion'
            ytry2=tryIt2(p,y,pave,mp,np,ndim,pstr,pstr2,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat,ihi,gam)
            iter=iter+1
c            print*,' tryIt2  ytry2 y(ilo) ', ytry2, y(ilo)
c   Box 5
 
            if (ytry2 .lt. y(ilo))  then
c               print*,' real improvement'
               y(ihi)=ytry2
               do 24 j=1,ndim
                 p(ihi,ifloat(j))=pstr2(ifloat(j))
24             continue
            else
c               print*,'  not improvement'
               y(ihi)=ytry
               do 25 j=1,ndim
                 p(ihi,ifloat(j))=pstr(ifloat(j))
25             continue
            end if
            itest=1
       end if
c  still part of box 1
       if (ytry .gt. y(ilo)) then
c  Box 2
         itest=0
         if ( (ytry .gt. y(ilo)) .and. (ytry .lt. y(ihi)) ) then
c            print*,' ytry between limits'
c   place yhi with ytry
            y(ihi)=ytry
            do 22 j=1,ndim
              p(ihi,ifloat(j))=pstr(ifloat(j))
22          continue
            itest=1
         else
            itest=0
c           print*, 'ytry is probably greater than y(ihi)'
c           print*,' ytry y_lo y_hi ', ytry,y(ilo),y(ihi)
c   Do a contraction
           ytry=tryIt(p,y,pave,mp,np,ndim,pstr,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat,ihi,beta)
c           print*,' contraction ytry and y_hi ', ytry,y(ihi)
           iter=iter+1
c  Box 4
           if (ytry .gt. y(ihi) )  then
           do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                pstr(ifloat(j))=0.5*(p(i,ifloat(j))+p(ilo,ifloat(j)))
                p(i,ifloat(j))=pstr(ifloat(j))
15            continue


              call funmin(YPRR,pstr(1),pstr(2),pstr(3),pstr(4),pstr(5),
     &    npole,flow,fhigh,pstr(6),pstr(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  
              iter=iter+1

            y(i)=-1.0*YPRR

            endif
16         continue
           itest=1
           else
             itest=1
             y(ihi)=ytry
             do 23 j=1,ndim
               p(ihi,ifloat(j))=pstr(ifloat(j))
23           continue

           end if
         end if

       end if
       if (itest .eq. 1) go to 1

      return
      END
      function tryIt(p,y,pave,mp,np,ndim,ptry,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat,ihi,fac)

      dimension d(max_data),res(max_data)
      double precision t_year(max_data), A(max_data,max_mod),dt_sam

      INTEGER ihi,mp,ndim,np,NMAX,ifloat(mp)
      REAL fac,p(mp,np),pave(np)
      PARAMETER (NMAX=20)

      INTEGER j
      REAL fac1,fac2,ptry(NMAX)
      double precision YPRR,ytry,tryIt,y(mp)
      fac1=(1.-fac)
c      fac2=fac1-fac

c      print*,"fac ",fac,fac1
c initializ ptry
      do 5 i=1,np
        ptry(i)=p(1,i)
5     continue
      do 11 j=1,ndim
        ptry(ifloat(j))=pave(ifloat(j))*fac1+p(ihi,ifloat(j))*fac
c        print*,j,ptry(j),pave(j),fac1,p(ihi,j),fac
11    continue
c1234567890123456789012345678901234567890123456789012345678901234567890123456789

c       print*,' 2',ptry(4)

        call funmin(YPRR,ptry(1),ptry(2),ptry(3),ptry(4),ptry(5),
     &    npole,flow,fhigh,ptry(6),ptry(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  


      ytry=-1.0d+0*YPRR

      tryIt=ytry
      return
      end
      function tryIt2(p,y,pave,mp,np,ndim,ptry,ptry2,
     &   npole,flow,fhigh,dt_sam,max_data,max_mod,
     &    t_year,iswt,modType,A,d,res,ifloat,ihi,fac)

      dimension d(max_data),res(max_data)
      double precision t_year(max_data), A(max_data,max_mod),dt_sam

      INTEGER ihi,mp,ndim,np,NMAX,ifloat(mp)
      REAL fac,p(mp,np),pave(np)
      PARAMETER (NMAX=20)

      INTEGER j
      REAL fac1,fac2,ptry(NMAX),ptry2(NMAX)
      double precision YPRR,ytry,tryIt2,y(mp)
      fac1=(1.-fac)
c      fac2=fac1-fac

c      print*,"fac tryIt2",fac,fac1
c initializ ptry
      do 5 i=1,np
        ptry2(i)=p(1,i)
5     continue
      do 11 j=1,ndim

        ptry2(ifloat(j))=pave(ifloat(j))*fac1+ptry(ifloat(j))*fac
c        print*,j,ptry(ifloat(j)),pave(j),fac1,fac,ptry2(ifloat(j))
11    continue
c1234567890123456789012345678901234567890123456789012345678901234567890123456789

c       print*,' 2',ptry(4)

        call funmin(YPRR,ptry2(1),ptry2(2),ptry2(3),ptry2(4),ptry2(5),
     &    npole,flow,fhigh,ptry2(6),ptry2(7),dt_sam/365.25d+0,
     &    iswt,ModType,A,d,res)  


      tryIt2=-1.0d+0*YPRR

      return
      END
