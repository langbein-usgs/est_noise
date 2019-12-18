      subroutine NedlerMeadEXP(p,y,MP,NP,NDIM,FTOL,ITER,
     &     ModType,A,d,res,max_data,max_mod)
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
c  Modified NedlerMead.f such that the function is the time-dependent model where
c    the time constant of the omori or exponential is being estimated
c
      dimension d(max_data),res(max_data)
c      double precision t_year(max_data),dt_sam,y(MP),yswap,tryItEXP2,
c     &   ybar,y2
      double precision y(MP),yswap,tryItEXP2,tryItEXP,
     &   ybar,y2
      double precision A(max_data,max_mod),tryIt2,ytry,YPRR,ysave,ytry2
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL ftol,p(mp,np), nave
      PARAMETER (NMAX=20)
      INTEGER i,ihi,ilo,inhi,j,m,n
      character*1 ModType
      REAL rtol,sum,swap,pave(NMAX),pstr(NMAX),pstr2(NMAX)

c  Note, my alpha is -1*alpha of nedler and mead --- allows single subroutine Pave*(1-val)+P*val to be evaluated
      alpha=-1.0
c If there is only one unknown, then nedler mead doesn't iterate... it becomes
c   degenerate; it automatically picks the best of two solutions,
c   so, to avoid this, make alpha not equal to -1.0
      if (ndim .eq. 1) alpha=-1.1
      beta=0.5
      gam=2.0
c  experiment
      alpha=-1.0 - 0.1/float(ndim)
      beta=0.5 + 0.05/float(ndim)
      gam=2.0+ 0.2/float(ndim)
      ITMAX=75
      iter=0
      rtolMin=999999.0
      nave=float(ndim+1)
c      print*,'ModType ',ModType
c      print*,' In NedlerMeadEXP, MP,NP,NDIM,FTOL,ITER',
c     *  MP,NP,NDIM,FTOL,ITER,
c     *  ' max_data max_mod', max_data,max_mod
c      do 200 n=1,ndim+1
c        print*,n,y(n),(p(n,i),i=1,ndim)
c200   continue
c      print*, ' above is imput to NedlarMead'
c      nmod=9
c      print*,' in NedlerMead output A and d, nmod=',nmod
c      print*,'d(100)', (99+k,d(99+k),k=1,10)
c      print*,(A(100,j),j=1,nmod)
c      print*,'d(100)',d(100)
1     do 12 n=1,ndim
        sum=0.
c        print*,(p(m,(n)),m=1,ndim+1)
        do 11 m=1,ndim+1
          sum=sum+p(m,(n))/nave
11      continue
        pave((n))=sum
c      write(89,*)' 1st paves',(n),pave((n))
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
c       write(89,*)' ihi ilo ',ihi,y(ihi),ilo,y(ilo)


c  measure deviations  and compare with ftol

c do absolute tolerance 
c      rtol=abs(y(ihi)-y(ilo))
c  extremes
      rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi)+y(ilo)))
c  RMS
c      ybar=0.0d+0
c      do 51 i=1,ndim+1
c        ybar=ybar+y(i)/nave
c51    continue
c      y2=0.0d+0
c      do 52 i=1,ndim+1
c        y2=y2+(y(i)-ybar)**2
c52    continue
c      rtol=(sqrt(y2)/nave)/ybar
      if (rtol .lt. rtolMin) rtolMin=rtol
c      print*,'rtol=',rtol,y(ihi),y(ilo),ihi,ilo
      if (rtol.lt.ftol) then
        yswap=y(1)
        y(1)=y(ilo)
        y(ilo)=yswap
        do 14 n=1,ndim
          swap=p(1,(n))
          p(1,(n))=p(ilo,(n))
          p(ilo,(n))=swap
14      continue
c        write(89,*)' Converged! iter= ',iter
        return
      endif
c      if (iter.ge.ITMAX) pause 'ITMAX exceeded in NedlerMeadeEXP'
       if (iter.ge.ITMAX) then
         PRINT*,  'NedlerMeadEXP exceeding maximum iterations.',ITMAX
         print*,' rtolMin= ',rtolMin,'rtol=',rtol,y(ihi),y(ilo)
         write(89,*)'NedlerMeadEXP exceeding maximum iterations.',ITMAX
         write(89,*)' rtolMin= ',rtolMin,'rtol=',rtol,y(ihi),y(ilo)
         return
       end if
      iter=iter+1
c    Do a reflection and compute y*
c        write(89,*)' Do reflection ytryItEXP'
        ytry=tryItEXP(p,y,pave,mp,np,ndim,pstr,
     &   ModType,A,d,res,ihi,alpha,max_data,max_mod)
c        write(89,*)' done with reflection  ytryItEXP', ytry

c   Box 1
       if (ytry .lt. y(ilo)) then
c   improvement
c         write(89,*)' this is an improvement tryItEXP2'
c      This is 'expansion'
            ytry2=tryItEXP2(p,y,pave,mp,np,ndim,pstr,pstr2,
     &        ModType,A,d,res,ihi,gam,max_data,max_mod)
            iter=iter+1
c            write(89,*)' iter ',iter
c            write(89,*)' tryItEXP2  ytry2 y(ilo) ', ytry2, y(ilo)
c   Box 5
 
            if (ytry2 .lt. y(ilo))  then
c               write(89,*)' real improvement'
               y(ihi)=ytry2
               do 24 j=1,ndim
                 p(ihi,(j))=pstr2((j))
24             continue
            else
c               write(89,*)'  not improvement'
               y(ihi)=ytry
               do 25 j=1,ndim
                 p(ihi,(j))=pstr((j))
25             continue
            end if
            itest=1
       end if
c  still part of box 1
       if (ytry .gt. y(ilo)) then
c  Box 2
         itest=0
         if ( (ytry .gt. y(ilo)) .and. (ytry .lt. y(ihi)) ) then
c            write(89,*)' ytry between limits'
c   place yhi with ytry
            y(ihi)=ytry
            do 22 j=1,ndim
              p(ihi,(j))=pstr((j))
22          continue
            itest=1
         else
            itest=0
c           write(89,*) 'ytry is probably greater than y(ihi)'
c           write(89,*)' ytry y_lo y_hi ', ytry,y(ilo),y(ihi)
c   Do a contraction
c           write(89,*), ' contraction tryItEXP2'
           ytry=tryItEXP(p,y,pave,mp,np,ndim,pstr,
     &      ModType,A,d,res,ihi,beta,max_data,max_mod)
c           write(89,*)' done with contraction tryItEXP2'
c           write(89,*)' contraction ytry and y_hi ', ytry,y(ihi)
           iter=iter+1
c            write(89,*)' iter ',iter
c  Box 4
           if (ytry .gt. y(ihi) )  then
           do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                pstr((j))=0.5*(p(i,(j))+p(ilo,(j)))
                p(i,(j))=pstr((j))
15            continue


              iter=iter+1

c            write(89,*)' start fmodel_fit_A'
            y(i)=fmodel_fit_A(0,ModType,A,d,res,pstr)
c            write(89,*)' done fmodel_fit_A'
c            write(89,*)'pstr =',(pstr(j),j=1,ndim)
c            write(89,*)'cycle through model, i,y(i)',i,y(i)
c            write(89,*)' iter ',iter
            endif
16         continue
           itest=1
           else
             itest=1
             y(ihi)=ytry
             do 23 j=1,ndim
               p(ihi,(j))=pstr((j))
23           continue

           end if
         end if

       end if
       if (itest .eq. 1) go to 1

      return
      END
      function tryItEXP(p,y,pave,mp,np,ndim,ptry,
     & ModType,A,d,res,ihi,fac,max_data,max_mod)


      dimension d(max_data),res(max_data)
      double precision t_year(max_data), A(max_data,max_mod),dt_sam

      INTEGER ihi,mp,ndim,np,NMAX
      REAL fac,p(mp,np),pave(np)
      PARAMETER (NMAX=20)
      character*1 ModType
      INTEGER j
      REAL fac1,fac2,ptry(NMAX)
      double precision YPRR,ytry,tryItEXP,y(mp)
      fac1=(1.-fac)
c      fac2=fac1-fac

c      write(89,*)"In tryItEXP fac ",fac,fac1
c initializ ptry
      do 5 i=1,np
        ptry(i)=p(1,i)
5     continue
      do 11 j=1,ndim
        ptry((j))=pave((j))*fac1+p(ihi,(j))*fac
c        write(89,*)j,ptry(j),pave(j),fac1,p(ihi,j),fac
11    continue
c1234567890123456789012345678901234567890123456789012345678901234567890123456789

c       print*,' 2',ptry(4)



c      write(89,*) ' In EXP start fmodel_fit_A'
      ytry=fmodel_fit_A(0,ModType,A,d,res,ptry)
c      write(89,*) ' In EXP done fmodel_fit_A'
      tryItEXP=ytry
      return
      end
      function tryItEXP2(p,y,pave,mp,np,ndim,ptry,ptry2,
     & ModType,A,d,res,ihi,fac,max_data,max_mod)

      dimension d(max_data),res(max_data)
      double precision  A(max_data,max_mod),dt_sam

      INTEGER ihi,mp,ndim,np,NMAX
      REAL fac,p(mp,np),pave(np)
      PARAMETER (NMAX=20)
      character*1 ModType
      INTEGER j
      REAL fac1,fac2,ptry(NMAX),ptry2(NMAX)
      double precision YPRR,ytry,tryItEXP2,y(mp)
      fac1=(1.-fac)
c      fac2=fac1-fac

c      write(89,*)"fac tryItEXP2",fac,fac1
c initializ ptry
      do 5 i=1,np
        ptry2(i)=p(1,i)
5     continue
      do 11 j=1,ndim

        ptry2((j))=pave((j))*fac1+ptry((j))*fac
c        write(89,*)j,ptry((j)),pave(j),fac1,fac,ptry2((j))
11    continue

c      write(89,*) ' In EXP2 start fmodel_fit_A'
      tryItEXP2=fmodel_fit_A(0,ModType,A,d,res,ptry2)
c      write(89,*) ' In EXP2 done fmodel_fit_A'
      return
      END
      function fmodel_fit_A(nopt,ModType,A,d,res,parm)
c max_data,max_parm,ic,nmod,A,t,d,res,
c     &  covinv,covar,iswitch,x,e,
c     &  nexp,bexp,texp,exp_choice,exp_type,t_start,t_stop,parm)

c  function called by Downhill simplex; the sum of square residuals are
c    output to AMOEBA_mod
c max_data  dimension of data
c max_parm  dimension of num parameters
c ic number of data
c nmod number of parameters
c ax  design matrix
c d  data
c t  time
c res  residuals
c covinv  inverse covariance matrix
c covar  covariance matrix
c iswitch  outputs the model and its error if set to 1
c  x  the model
c  e  the error
c  nexp  number of exponentials both "fix" and "float"
c  bexp  are the time constants for the exponential (floaters are trial time constrants
c  texp  are the T_o of each of the exponentials
c  exp_choice is either "fix" or "float"
c  t_start  time, in days of first pt in time series
c  t_stop    time, in day of last pt in time series
c   parm is the trial time constants for expontials/omori
c
      double precision A(max_data,max_mod),dtime(11687),
     & covinv(11687,11687),covar(11687,11687),chi2
      double precision t_start,t_stop,texp(82),sumRes,t_year(11687)
      dimension d(max_data),res(max_data),x(82),e(82),
     & bexp(20),parm(20)
      character*1 ModType
      character*7 exp_choice(10),exp_type(10)
c      dimension sres(11),spar(11,10)

c      integer irowmiss(11687)

c   double precision
      common /ModFit1/dtime,t_year,covinv,covar,texp
      common /ModFit1a/t_start,t_stop
c  single precision
      common /ModFit2/x,e,bexp
c integer
      common /ModFit3/max_data,max_mod,ic,nmod,n_exp,nmiss,irowmiss,
     & max_time,ipl_flag_1,ipl_flag_2,ibp_flag,nmiss_max
c  character
      common /ModFit4/ exp_choice,exp_type
      mdmax=max_data
      maxnmod=max_parm
      nmod_orig=nmod
      nexp=n_exp
c      print*,' In function fmodel_fit'
c      print*," nexp=",nexp," ic=",ic," t_start=",t_start
c      print*,(k,exp_choice(k),k=1,nexp)
c  For the initial solution
      ix=0
c      print*," "
c      print*," "
      do 1 i=1,nexp
c      print*,i,texp(i),bexp(i)
      if (exp_choice(i) .eq. "float") then
         ix=ix+1
c         if (parm(ix) .lt. -4.) parm(ix)=-4.0
c         if (parm(ix) .gt. 2.0) parm(ix)=2.0
         bexp(i)=10.0**parm(ix)

c         print*, "Try ",ix,parm(ix),bexp(i)
      end if
1     continue
c      print*,' in fmodel_fit output A and d'
c      print*,(A(100,j),j=1,nmod)
c      print*,d(100),dtime(100)
      call modify_A(A,dtime,texp)
c      call calcres(0,nopt,sumRes,A,d,res,chi2)
      call calcres(nopt,ModType,sumRes,A,d,res,chi2)
c      write(13,130)(x(i),e(i),i=1,nmod)
c130   format(30(3x,f15.3,1x,f15.3))
c      print*, "finish calc_res",(x(i),i=1,nmod)


c      print*, " Chi^2 is", sum," RMS=",sqrt(sum/float(ic-nmod))
      nmod=nmod_orig
      fmodel_fit_A=chi2

      write(89,*)chi2,(parm(ix),10.0**parm(i),i=1,ix)
        
      return
      end
c
c


