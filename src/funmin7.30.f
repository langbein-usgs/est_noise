      subroutine funmin(fmle,sig1,sig2,sig3,sig4,sig5,
     &   npole,flow,fhigh,sig6,sig7,ts_yrs,
     &   iswt,modType,ax,d,res)

c
c  computes maximum likelihood for specified inputs
c
c  fmle Log maximum likelihod
c  sig1  white-noise
c  sig2  random walk noise
c  sig3  exponent of PSD
c  sig4  time-constant of PSD
c  sig5 is amplitude of Band-pass filtered noise
c    npole is the number of poles for BP filter
c    flow and fhigh is the pass-band in c/yr for BP filter
c  sig6  amplitude of second power law noise
c  sig7  exponent of second power law
c  ts_yrs  sampling interval in years
c  t_year  time in days of data (double precision
c  d   data is sequential order
c  txx  time of data obs in years
c  jmax  number of observations
c  max_data and max_parm  dimensions of matrices
c  ax   design matrix
c  nmod  number of parameters in design matrix
c  ipl_flag_1  a flag; if 1, power law covariance matrix has been computed for
c  iswt  switch for yes/no calc of error bars for model
c  irowmiss  --- indices where there are missing data
c  nmiss -  number of missing data;  note that the A and d matrix/vector are zero for each missing term
c
c      call funmin(fmle,siginstr,amp1,exp1,alpha,amp_bp,
c     &   npole,flow,fhigh,amp2,exp2,dt_sam/365.25,
c     &  iswt)
c  April 2021 --- modify the GGM model to use cholesky decomposition when iswt is 1 (output lots of stuff)
c    Does this by inserting 'missing' data prior to the start of the time-series to
c     create the proper covariance matrix for GGM model

      dimension d(max_data),x(82),xe(82),
     & res(max_data),
c     &  Fm(12687,12687),
     &  dxx(max_data),corrxx(max_data)
      dimension dr(12687)
      double precision dsum,t_year(12687),
c     & FtCi(12687,12687),FtCiF(12687,12687),
     & AAA(nmiss_max,nmiss_max),axx(max_data,max_parm),
     & AAAA(max_data,nmiss_max),ax(max_data,max_parm),chi2,fmle

      double precision dett,detmiss,covar(12687,12687),
     &  covinv(12687,12687)
      double precision filtwn(32768),filtpl1(32768),
     & filtpl2(32768),filtbp(32768)
      double precision covarpl1(12687,12687),covarpl2(12687,12687),
     & covarbp(12687,12687)

      double precision f(max_time),finv(32768), ts_yrs,small,fsmall,
     & t_small
      double precision B(12687,12687),E(max_data,nmiss_max)

      double precision texp(10),t_start,txx(12687),t_stop
      dimension bexp(10), parm(10)
      character*7 exp_choice(10), exp_type(10)
      integer irow(12687), irowmiss(12687)
      character*1 modType
c  Next used for DFT
c      double precision ftemp(65536),wsave(196623)
c   double precision
      common /ModFit1/txx,t_year,covinv,covar,texp
      common /ModFit1a/t_start,t_stop
c  single precision
      common /ModFit2/x,xe,bexp
      common /ModFit2a/expmax
c integer
      common /ModFit3/max_data,max_parm,jmax,nmod,n_exp,nmiss,irowmiss,
     & max_time,ipl_flag_1,ipl_flag_2,ibp_flag,nmiss_max
c  character
      common /ModFit4/ exp_choice,exp_type
c  Intermediate covariance and filter functions; double precision
      common /CovFlt1/filtpl1,filtpl2,filtbp,covarpl1,covarpl2,covarbp
c  double precision
      common /Modfit5/ finv, B

      md=max_data
      maxnmod=max_parm
      mmax=max_time
      irowOffset=0
      nexp=n_exp
c      print*,'Modtype=', Modtype

      tfunf0=secnds(0.0)
c      print*,' start funmin'
c       print*,jmax,max_data,max_parm,max_time
c       print*,fmle,sig1,sig2,sig3,sig4,sig5,
c     &   npole,flow,fhigh,sig6,sig7,ts_yrs
c       print*,' time variables',jmax,
c     & t_year(jmax),t_year(1),txx(jmax),txx(1)
c       print*,' nmiss=',nmiss, t_start,t_stop
c       do 998 i=1,jmax
c998    write(99,*),i,t_year(i),txx(i)
c       print*,covinv(1,1),covar(1,1),texp(1),texp(10),t_start,t_stop
 

c      print*,'funmin',sig1,sig2,sig3,sig4,sig5,sig6,sig7
      iswitch_verbose=0


c   if iswitch_verbose=1----lots of output

3890  continue
      if (jmax .gt. md) then
       print*, 'Number of observations, ',jmax,' greater than',
     &         ' storage ',md
       stop
      end if
      if (nmod .gt. maxnmod) then
       print*, 'Number of model param., ',nmod,' greater than',
     &   ' storage ',maxnmod
       stop
      end if
      small=1.0d-15
c  make sure that everything is positive except  but exponent to power law
      sig1=abs(sig1)
      sig2=abs(sig2)
      sig4=abs(sig4)
      sig5=abs(sig5)
      sig6=abs(sig6)
      sig7=abs(sig7)
      if (sig2 .lt. 0.0001) sig2=0.0001
      if (sig6 .lt. 0.0001) sig6=0.0001
c  limit the first index to expmax
      if (sig3 .gt. expmax) sig3=expmax
c  
c
c  figure-out the equivalent PSD parameters....the Po term!
c

      fs=1./ts_yrs
      t_small=ts_yrs
      if (sig2 .ne. 0 ) then
c        pow_db=5*sig3*alog10(fs)-7.339*sig3-10*alog10(fs)+2.076
        Po=2.0*(6.28**(-1.0*sig3))/(fs**(1-(sig3/2)))*sig2*sig2
c        pow_db=20*alog10(sig2)+pow_db
c        Po=10**(psd_db/10.)
        pow_db=10.0*alog(Po)
      end if

c
c  compute ratio to relate Po to the random walk number
c
c
c  figure-out lowest value of sig1 detectable
c
c  compute smallest sampling interval

      fny=1/(2*t_small)
      sig1_orig=sig1
      fs=1./ts_yrs
      sig_min=sig2*sqrt(2.0)*((2./3.14159)**(sig3/2.))/(fs**(sig3/4))

c      print*," sig1=", sig1, " sig_min=", sig_min
c      if (sig1 .lt. sig_min/2.) sig1=sig_min/2.0
c      print*, " sig1=", sig1


      if (sig4 .lt. 0) then
        print*,' '
        print*,' Time constant was ',sig4
        sig4=0.0001
        sig4=2.0*3.14159/(60*t_year(jmax))
        print*,' force time constant to be ',sig4
      end if
c      if (sig4 .gt. 6.28/t_year(jmax))
c     &     sig4=6.28/t_year(jmax)
      sig4_orig=sig4
      

      siginstr=sig1**2
c test for negative amplitude for band-pass filter noise
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
        print*,' amplitdue of PSD at lowest frequency (DC) in db is ',
     &  psd_db -10*sig3*alog10(fo)
      end if
      if (sig1 .ne. 0) then
c        print*,fny,sig1
        print*,' amplitude of PSD at highest frequency is db is: ',
     &  10*alog10(sig1*sig1/fny)
        fc=(((Po)*fny)/(sig1*sig1))**(1./sig3)
        print*,' cross over freq (c/yr) ',fc
      end if
      print*,' BP-filter amplitude: ',sig5
      print*,' BP-filter parameters; Poles= ',npole,
     & ' low and high freqs; ',flow,fhigh,' c/yr'
      if (sig6 .ne. 6) then
        print*,' Amplidude of second power law noise ',sig6,
     &   ' units/yr^0.25*sig7'
        print*,' Exponent of second power law noise', sig7

      end if

      end if

       time0=secnds(0.0)


       kmax=int((t_year(jmax)-t_year(1))/t_small + 0.5)+1
       if (kmax .gt. max_time) then
         print*,' kmax = ',kmax,' exceeds max_time of ',max_time
         stop
       end if

      if ((ModType .eq. 'f' ) .or. (ModType .eq. 'c')) then

c  create first power law and GM filter

      if (ipl_flag_1 .ne. 1)  then
        kkmax=kmax
        irowOffset=0
c        if ((iswt .eq. 1) .and. (sig4 .gt. 0.0)) then
        if ( (sig4 .gt. 0.0) .and. (ModType .eq. 'c')) then
c  normally, create 'ghost' data preceding actual time series for GM
c     for GM noise
c   But, for sig4 roughly equal to the length of the time series, it
c    is close to PL noise and those ghost data aren't needed to compute
c    an accurate covariance matrix
          ttlen=(t_year(jmax)-t_year(1))
          fxlow=1/ttlen
          extend=3.0
c  compare sig4 to ttlen
         fnaught=sig4/6.28
         if (fnaught .lt. fxlow) then
           extend=-1.5*(fxlow-fnaught)+3
           if (fnaught .lt. 0.5*fxlow) extend=0
         end if
c   check to see if sig4 could make covariance matrix exceed dimensions
          if (sig4 .lt. 3/((float(max_time-kmax)*sngl(t_small)) ))
     &       sig4=3.10/((float(max_time-kmax)*sngl(t_small)) )
c  extend time series by 3/sig4 (sig4 is ggm in rad/year)
c           print*,' orig kkmax=', kkmax,t_small,sig4
           kkmax=kmax+ifix(extend*1/(sngl(t_small)*sig4))
c           print*,sig4,kkmax,kkmax-kmax,kmax, max_time
c           print*,fxlow,fnaught," extend", extend
c  what happens when sig4 is very close to zero and requests something
c    that exceeds array size..... put limit on kkmax.
           if (kkmax .gt. max_time) kkmax=max_time
           irowOffset=kkmax-kmax
c            print*,' kkmax=', kkmax, 'irowOffset=',irowOffset
           if (kkmax .gt. max_time ) then
              print*,' kmax = ',kmax,' exceeds max_time of ',max_time
              stop
           end if           
        end if
        call frac_diff(filtpl1,sig3,sig4,sngl(t_small),kkmax,mmax)
      end if
      const1=(sig2**2)*(t_small**(sig3/2.0)) 
      const1=sqrt(const1)

c  create second power law
      if (ipl_flag_2 .ne. 1) then
        call frac_diff(filtpl2,sig7,0.0,sngl(t_small),kmax,mmax)
      end if
      const2=(sig6**2)*(t_small**(sig7/2.0)) 
      const2=sqrt(const2) 

c  Note --- is ibp_flag is always 1, then this is useless code

       if (ibp_flag .ne. 1) then
       if (sig5 .ne. 0) 
     &  call band_pass_cov(t_year,jmax,t_small,flow,fhigh,npole,
     &   mmax,md,1.0,covarbp)
c     &   print*,' do bp'
       end if


c
c Add all of the filters
c
       do 10 i=1,kmax
         f(i)=const1*filtpl1(i) + const2*filtpl2(i)+ sig5*filtbp(i)

10     continue

       f(1)=f(1)+sig1
       if ( ModType .eq. 'f' ) then
c  DFT approach to de-convolving the filter --  Commented OUT (including a second end if
c       if (sig5 .ne. 0 ) then
c         print*,' revise filter function by computed PSD and smoothing'
c   time shift data and cosine taper
c         pi=3.141592653589793
c         do 120 i=1,kmax
c           ftemp(i)=0.0
c           taper=cos(2*pi*float(i)/(4.0*float(kmax)))
cc           taper=taper**2
cc           taper=1.0
c           ftemp(i+kmax)=f(i)*taper
c           f(i)=f(i)*taper
cc           write(44,*)i,ftemp(i),ftemp(i+kmax),taper
c120      continue
c         knn=kmax*2
        
c         call dffti(knn,wsave)  
c         call dfftf(knn,ftemp,wsave)

c         do 121 i=2,kmax
c           dre=ftemp(i*2-2)/float(knn)
c           ddi=ftemp(i*2-1)/float(knn)
c           ftemp(i*2-2)=dre/((dre**2 + ddi**2)*(float(knn))**2)
c           ftemp(i*2-1)=-1.0*ddi/((dre**2 + ddi**2)*(float(knn))**2)    
            
c121      continue


c         ftemp(1)=ftemp(1)/float(knn)
c         ftemp(1)=ftemp(1)/(ftemp(1)**2*float(knn)**2)
c         ftemp(knn)=ftemp(knn)/float(knn)
c         ftemp(knn)=ftemp(knn)/(ftemp(knn)**2*float(knn)**2)
c         call dfftb(knn,ftemp,wsave)

c         do 123 i=1,kmax
cc           write(41,*)i,ftemp(i),ftemp(i+kmax)
cc           finv(i)=ftemp(i+kmax)-ftemp(i)
c           finv(i)=ftemp(i+kmax)
          
c123      continue
c       else

c
c  invert the filter
c    by inverse convolution
      fsmall=1.0d-30
      finv(1)=1.0/f(1)
c      write(45,*)1,f(1),finv(1)
      do 15 i=2,kmax
        dsum=0.d+0
        do 14 j=1,i-1
          dsum=dsum+(finv(j))*(f(i+1-j))
14      continue
        finv(i)=-1.0*(dsum/(f(1)))
c  Needed --- if finv is zero, the construction of covariance matrix takes too long;
        if (abs(finv(i)) .lt. fsmall) finv(i)=fsmall
c      write(45,*)i,f(i),finv(i)


15    continue

C commented out 1st end if for DFT
c      end if



c
c Test filter
c
      iiverb=1
      if (iiverb .eq. 0) then
      dsum=0.0d+0
      do 18 i=1,kmax
        fmle=0.0d+0
        do 17 j=1,i
          fmle=fmle+f(j)*finv(i-j+1)
c          fmle=fmle+finv(j)*f(i-j+1)
17      continue
        if (i .ne. 1) dsum=fmle**2+dsum
        if (i .eq. 1) dsum=(1.0d+0-fmle)**2 + dsum
c        write(46,*)i,fmle

18    continue
      print*,' Inverse test',sqrt(dsum/(kmax))
      end if
c      if (kmax .eq. kmax) stop
c
c  get indices when f(i) and finv(i) are small
c
      ifflag=0
      ifinvflag=0
      do 181 i=jmax,1,-1
        if (ifflag .eq. 0) then
          if (abs(f(i)) .gt. fsmall) then
            ifflag=1
            maxf=i
          end if
        end if
        if (ifinvflag .eq.0) then
          if (abs(finv(i)) .gt. fsmall) then
            ifinvflag=1
            maxfinv=i
          end if
       end if
181   continue
c  end for Modtype .eq. 'f'
      end if

c   if ModType=f  use bos et al
c   if ModType=c  do cholesky decomposition of covariance (ala version 7.02)
      if (ModType .eq. 'c') then

c  Do cholesky to get inverse

c
c   Revise the data and design matrix by tossing out rows of missing data
c
       ix=1
       ixx=1
       do 8001 i=1,nmiss
         krow=irowmiss(i)
         do 8002 k=ix,krow-1
           dxx(ixx)=d(k)
           corrxx(ixx)=t_year(k)
c           print*,'i,krow,ixx,k',i,krow,ixx,k
           do 8003 j=1,nmod
             axx(ixx,j)=ax(k,j)
8003       continue           
           ixx=ixx+1
8002    continue
c        print*,"krow,d(krow),t_year(krow)",krow,d(krow),t_year(krow)
        ix=krow+1

8001   continue
       do 8004 k=ix,jmax
         dxx(ixx)=d(k)
         corrxx(ixx)=t_year(k)
c         print*,'ixx,k',ixx,k
         do 8005 j=1,nmod
           axx(ixx,j)=ax(k,j)
8005     continue
         ixx=ixx+1
8004   continue
c       print*,'jmax,kmax,nmiss,ixx',jmax,kmax,nmiss,ixx
       ixx=ixx-1
       do 2011 i=1,ixx

           ix=int((corrxx(i)-corrxx(1))/t_small+0.5)+1
           irow(i)=ix+irowOffset
c           print*,i,ix

2011  continue


      call makeMatrix(irow,kmax+irowOffset,f,covar,ixx,0,mmax,md)

      do 2801 i=1,ixx
        do 2802 j=1,i
        covinv(i,j)=covar(i,j)
        covar(j,i)=covar(i,j)
        covinv(j,i)=covar(i,j)

2802    continue

2801  continue




       call dpotrf('L',ixx,covinv,md,nerror)

c  bail-out if nerror not equal 0
      if (nerror .ne. 0) then
      print*,' covariance is singular, nerror=',nerror
c      print*,sig1,sig2,sig3,sig4,sig5,sig6,sig7

      
      sig1_new=sig1*1.26
      if (sig1 .lt. 0.25) sig1_new=0.40
      print*,' replace amplitude of noise,', sig1,' by ',sig1_new
      sig1=sig1_new
      if (sig4 .ne. sig4_orig) sig4=sig4_orig
      if (sig2 .ne. sig2_orig) sig2=sig2_orig
c      if (nerror .ne. 0) stop
      go to 3890
      end if
c  copy lower triangular part of covinv
c  calc the determanant
c  search for stuff close to zero

      dett=0.0
      covar_max=0.
      covar_min=1.0e+20
      covar_minx=1.0e+20
      do 730 i=1,ixx
        dsum=covinv(i,i)
        dett=dett+dlog(dsum)
        do 7328 j=1,i
        if (abs(covinv(i,j)) .gt. covar_max) then
          covar_max=covinv(i,j)
          i_max=i
          j_max=j
        end if
        if (abs(covinv(i,j)) .lt. covar_min) then
          covar_min=covinv(i,j)
          i_min=i
          j_min=j
        end if 
        if (  (abs(covinv(i,j)) .lt. covar_minx ) .and.
     &       (abs(covinv(i,j)) .ne. 0.0 ) ) then
          covar_minx=covinv(i,j)
          i_minx=i
          j_minx=j
       end if 
7328   continue  
730   continue
       c_large=1.0e+30
       if (covar_max/covar_minx .gt. c_large) then
         covmin=covar_max/c_large
         n0=0
         do 70063 i=1,ixx
         do 70064 j=1,i
c            if (covinv(i,j) .ne. 0.0) then
            if (abs(covinv(i,j)) .lt. covmin) then
              covinv(i,j)=covmin

              n0=n0+1
c            end if
            end if
70064    continue
70063    continue
       end if
       call dpotri('L',ixx,covinv,md,nerror)


c      print*,' spotri, nerror',nerror, jmax,md
c  fill in upper trianguangular part of matrix
c
      do 7325 i=1,ixx

      do 7324 j=i+1,ixx
      covinv(i,j)=covinv(j,i)

7324  continue
7325  continue
      nmissx=0

      jmaxTmp=jmax
      jmax=ixx

      call model_fit(iswt,ModType,axx,dxx,dr,chi2)
c      call model_fit(max_data,max_parm,ixx,nmod,axx,txx,dxx,dr,
c     &  covinv,covar,iswt,x,xe,
c     &  nexp,bexp,texp,exp_choice,exp_type,t_start,t_stop,
c     &  nmissx,irowmiss)
      jmax=jmaxTmp
      end if


C ___________Above --- compute covariance old way using cholesky decomposition on large matrix
C -----------Below, use Bos et al. -------
      if (ModType .eq. 'f')  then

c      tt0=secnds(0.)
c
c  construct eqn 7 of Bos et al
c   C^-1 - C^-1 F ( F^t C^-1 F )^-1 F^t C^-1
c   where F represents the covainace of missing data (the columns in C for which data are missing)
c    C is covariance matrix for which there is no missing data
c
c  Revise the above;
c   Let F be the matrix composed of finv (inverse filter)
c   Let M be the matrix composed of indices for missing data
c   C^-1 = F^t F
c   Then, let E=F M
c  Re-write eqn 7 of Bos et al interms of new variables;
c   C^-1 - C^-1 M ( M^t C^-1 M )^-1 M^t C^-1
c  plug in C^-1 = F^t F
c  F^t F - F^t F M [ M^t F^t F M ]^-1 M^t F^t F
c  F^t F - F^t E [ E^t E ]^-1 E^t F

c  Let B = E [ E^t E ]^-1 E^t
c then
c  F^t F - F^t B F
c
c  When solving for model, x;
c  x = [A^t C^-1 A ] A^t C^-1 d
c
c  let Aw= F A
c      dw= F d   Whitening of both model and data;
c
c  x = [ Aw^t Aw - Aw^t B Aw]^-1 [ Aw^t - Aw^t B ] dw
c
c  So, in this portion, I will calculate B and the determinant;
c   which is n*log(f(1) + trace [ E^t E ]^-1

       covmax=0.0
       covmin=1.0e+20
      do 201 i=1,kmax
        irow(i)=i
201   continue



c
c   the determinant  --- 
       dett=float(jmax)*dlog(f(1))


c      if ( iswt .eq. 1) then
c   only needed to compute standard errors to the model
c        print*,' create covar'
 
c         call makeMatrix(irow,kmax,f,covar,kmax,0,mmax)



c        call makeMatrix(irow,kmax,finv,covinv,kmax,2,mmax)
c      end if

c      if (iswt .ne. 1) then
      if (iswt .ne. 2) then

      if (nmiss .ne. 0 ) then

c  Construct the E matrix; E=F M  where F is the inverse filter matrix and M are columns showing missing data indices
c

      do 2071 i=1,nmiss
        do 208 j=1,irowmiss(i)-1
        E(j,i)=0.0

208     continue
        do 2081 j=1,kmax-irowmiss(i)+1
        k=irowmiss(i)+j-1
        E(k,i)=finv(j)

2081    continue
2071  continue

c
c  compute E^t E
c

      call dgemm('T','N',nmiss,nmiss,kmax,1.0d+0,E,md,E,md,
     &  0.0d+0,AAA,nmiss_max)




c  Invert AAA=E^t E using cholesky routines
      call dpotrf('L',nmiss,AAA,nmiss_max,nerror)


c  Compute determinant due to missing data --- works
c   in conjunction with section below that does cholesky decomposition
c    of inverse covariance matrix with time series having no gaps...
c    This is replaced by the approximation above that only uses f(1) to get determinant
      detmiss=0.0d+0
      do 2092 i=1,nmiss
        detmiss=detmiss+dlog(AAA(i,i))
2092  continue
c      print*,' determinant for missing data', detz

      dett=dett+detmiss

c
c   continue with inversion
       call dpotri('L',nmiss,AAA,nmiss_max,nerror)

c  fill in upper trianguangular part of matrix
      do 207 i=1,nmiss
        do 206 j=i+1,nmiss
        AAA(i,j)=AAA(j,i)
206     continue
c        write(66,5555)(AAA(i,k),k=1,nmiss)
207   continue


c     
c    compute C^-1*F*(F^t*C^-1*F)^-1
c  compute B = E [ E^t E ]^-1 E^t = E AAA E^t
c
c    E AAA = AAAA

c      call dgemm('N','N',kmax,nmiss,kmax,1.0d+0,E,md,AAA,md,
c     &  0.0d+0,AAAA,md)
      call dsymm('R','L',kmax,nmiss,1.0d+0,AAA,nmiss_max,E,md,
     & 0.0d+0,AAAA,md)

c   B = E AAA E^t = AAAA E^t

      call dgemm('N','T',kmax,kmax,nmiss,1.0d+0,AAAA,md,E,md,
     &  0.0d+0,B,md)


c  end of if statement to test whether there are missing data
      end if

      else
c
c  compute inverse covariance the slower way --- used for computing error bars
c   pick rows of C^-1
c
       print*, ' iswt = ', iswt
ccc comment out this section!!  use ccc
ccc       do  3071 i=1,nmiss
ccc         do 308 j=1,kmax
ccc         FtCi(i,j)=covinv(irowmiss(i),j)
ccc308      continue
ccc3071   continue 


C  F^t*C^-1*F = FtCi*F = FtCiF


c  pick columns of F^t*C^-1 matrix
ccc      do 209 i=1,nmiss
ccc        do 2091 j=1,nmiss
ccc        FtCiF(i,j)=FtCi(i,irowmiss(j))

ccc2091    continue

ccc209   continue


c  Invert FtCiF using cholesky routines
ccc      call dpotrf('L',nmiss,FtCiF,md,nerror)


c  Compute determinant due to missing data --- works
c   in conjunction with section below that does cholesky decomposition
c    of inverse covariance matrix with time series having no gaps...
c    This is replaced by the approximation above that only uses f(1) to get determinant
ccc      detmiss=0.0d+0
ccc      do 3092 i=1,nmiss
ccc        detmiss=detmiss+dlog(FtCiF(i,i))
ccc3092  continue

ccc      dett=dett+detmiss
c
c   continue with inversion
ccc       call dpotri('L',nmiss,FtCiF,md,nerror)
c     
c    compute C^-1*F*(F^t*C^-1*F)^-1

c   (F^t*C^-1*F)^-1*(F^t*C^-1)  
ccc       call dsymm('L','L',nmiss,kmax,1.0d+0,FtCiF,md,FtCi,md,
ccc     &   0.d+0,AAA,md)

c
c    compute C^-1*F*(F^t*C^-1*F)^-1*F^t*C^-1

c  and simultaneously subtract from C^-1



ccc       call dgemm('T','N',kmax,kmax,nmiss,-1.0d+0,FtCi,md,AAA,md,
ccc     & 1.0d+0,covinv,md)

c  end ccc comment-out

      end if
c  end statement for iswt



        call model_fit(iswt,ModType,ax,d,dr,chi2)


      end if
C  ------Above --- Bos et al covariance




c      timeZ=secnds(0.0)
c      print*,' model done',timeZ-time0  
c
c multiply inv(chol(covinv))*d and square the cumulative 
c
c  need to do transpose!

c      sumres=0.
c      if (ModType .eq. 'f' ) then
cc  do this for Bos et al Covariance
c      do 7320 i=1,jmax
c      res(i)=0.
c          do 7321 j=1,jmax
c          res(i)=res(i)+dr(i)*covinv(i,j)*dr(j)
c7321      continue
c      sumres=sumres+res(i)
c7320  continue
c      end if
c      if (ModType .eq. 'c' ) then
c  do this for cholesky version
c      do 73201 i=1,jmax
c      res(i)=0.
c          do 73211 j=1,jmax
c          res(i)=res(i)+dr(i)*covinv(i,j)*dr(j)
c73211      continue
c      sumres=sumres+res(i)
c73201 continue
c      end if

      end if

c   End Modtype=f or c

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c  Do this if ModType is 'n'  OR STANDARD COVARIANCE!
c
       if (iswitch_verbose .eq. 1) then
       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' start computing covariances', del_sec
       end if
      if (ModType .eq. 'n' ) then
      if (ipl_flag_1 .ne. 1)  then

         call pow_law_cov(t_year,jmax,sig1,sig2,sig3,sig4,
     &    covar,mmax,md,t_small)

      end if
       if (iswitch_verbose .eq. 1) then
       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' first covariances', del_sec
       end if
      if (ipl_flag_2 .ne. 1) 
     &  call pow_law_cov(t_year,jmax,0.0,1.0,sig7,0.0,
     &  covarpl2,mmax,md,t_small)
c
c construct covariance do to BP-filtered noise
c
       if (ibp_flag .ne. 1) then
       if (sig5 .ne. 0) 
     &  call band_pass_cov(t_year,jmax,t_small,flow,fhigh,npole,
     &   mmax,md,1.0,covarbp)
       end if
       if (iswitch_verbose .eq. 1) then
       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' start adding covariances', del_sec
       end if
       if ((ipl_flag_1 .eq. 1) .and. (sig2 .ne. 0)) then
       do 7000 i=1,jmax
       do 7001 j=1,i-1
       covar(i,j)=(sig2**2)*covarpl1(i,j)
7001   continue
       covar(i,i)=(sig2**2)*covarpl1(i,i)+sig1**2
7000   continue
       end if
       if (sig5 .ne. 0) then
       do 7002 i=1,jmax
       do 7003 j=1,i
       covar(i,j)=(sig5**2)*covarbp(i,j)+covar(i,j)
7003   continue
7002   continue
       end if

       if (sig6 .ne. 0) then
       do 7004 i=1,jmax
       do 7005 j=1,i
       covar(i,j)=(sig6**2)*covarpl2(i,j)+covar(i,j)
7005   continue
7004   continue
       end if

       covar_max=0.
       covar_min=1.0e+20
       covar_minx=1.0e+20
       do 7006 i=1,jmax

       do 7007 j=1,i-1
       if (abs(covar(i,j)) .gt. covar_max) then
          covar_max=covar(i,j)
          i_max=i
          j_max=j
       end if
       if (abs(covar(i,j)) .lt. covar_min) then
          covar_min=covar(i,j)
          i_min=i
          j_min=j
       end if 
       if (  (abs(covar(i,j)) .lt. covar_minx ) .and.
     &       (abs(covar(i,j)) .ne. 0.0 ) ) then
          covar_minx=covar(i,j)
          i_minx=i
          j_minx=j
       end if   
       covar(j,i)=covar(i,j)
       covinv(i,j)=covar(i,j)
       covinv(j,i)=covar(i,j)
7007   continue
       covinv(i,i)=covar(i,i)
c       write(99,*)i,covinv(i,i)
7006   continue
       if (iswitch_verbose .eq. 1) then
         print*,' max covar = ',covar_max,' at element',i_max,j_max
         print*,' min covar = ',covar_min,' at element',i_min,j_min
         print*,' min covar0 = ',covar_minx,' at element',i_minx,j_minx
       end if

      c_large=1.0e+30

c
c  do cholsky
c

       del_sec=timex-time0
       
       if (iswitch_verbose .eq. 1) then

       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' start spotrf', del_sec
       end if
       nerror=0
c       print*,'jmax, md',jmax,md
       call dpotrf('L',jmax,covinv,md,nerror)
c       print*,' dpotrf nerror=',nerror
c       do 6003 i=1,5
c6003   print*,(covinv(i,j),j=1,5)
c       do 6004 i=jmax-5,jmax
c6004   print*,(covinv(i,j),j=jmax-5,jmax)
c       print*,(i,covinv(i,i),i=1,jmax)

       if (iswitch_verbose .eq. 1) then

       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' finish spotrf', del_sec
       end if
c  bail-out if nerror not equal 0
      if (nerror .ne. 0) then
      print*,' covariance is singular'
      sig2_new=sig2*0.78
      print*,' replace amplitude of noise,', sig2,' by ',sig2_new
      sig2=sig2_new
      if (sig4 .ne. sig4_orig) sig4=sig4_orig
      if (sig1 .ne. sig1_orig) sig1=sig1_orig
      go to 3890
      end if
      
c  copy lower triangular part of covinv
c  calc the determanant
c  search for stuff close to zero

      dett=0.0
      covar_max=0.
      covar_min=1.0e+20
      covar_minx=1.0e+20
      do 7300 i=1,jmax
        dsum=covinv(i,i)
        dett=dett+dlog(dsum)
c        write(99,*),i,covinv(i,i),dett
c        do 7329 k=1,jmax
c7329   covinv(i,k)=covinv(i,k)
        do 73280 j=1,i
        if (abs(covinv(i,j)) .gt. covar_max) then
          covar_max=covinv(i,j)
          i_max=i
          j_max=j
        end if
        if (abs(covinv(i,j)) .lt. covar_min) then
          covar_min=covinv(i,j)
          i_min=i
          j_min=j
        end if 
        if (  (abs(covinv(i,j)) .lt. covar_minx ) .and.
     &       (abs(covinv(i,j)) .ne. 0.0 ) ) then
          covar_minx=covinv(i,j)
          i_minx=i
          j_minx=j
       end if 
73280  continue  
7300  continue
      if (iswitch_verbose .eq. 1) then
        print*, 'Natural log of determinant ', dett
         print*,' max covar = ',covar_max,' at element',i_max,j_max
         print*,' min covar = ',covar_min,' at element',i_min,j_min
         print*,' min covar0 = ',covar_minx,' at element',i_minx,j_minx
 
         timex=secnds(0.0)
         del_sec=timex-time0
         print*,' start spotri',del_sec
      end if
c      c_large=1.0e+30
       if (covar_max/covar_minx .gt. c_large) then
         covmin=covar_max/c_large
         n0=0
         do 71063 i=1,jmax
         do 71064 j=1,i
c            if (covinv(i,j) .ne. 0.0) then
            if (abs(covinv(i,j)) .lt. covmin) then
              covinv(i,j)=covmin

              n0=n0+1
c            end if
            end if
71064   continue
71063   continue
           if (iswitch_verbose .eq. 1) then
             print*,' Number of elements replace by 0 is',n0
             timex=secnds(0.0)
             del_sec=timex-time0
             print*,' Zeroed-out covariance', del_sec
           end if
         end if
c
c invert cholesky matrix
c

       call dpotri('L',jmax,covinv,md,nerror)

       if (nerror .ne. 0) then
          print*, ' finished spotri but result is singular'
          print*,'   nerror=',nerrror
          stop
       end if
       if (iswitch_verbose .eq. 1) then

       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' finish spotri', del_sec
       end if
c  fill in upper trianguangular part of matrix
c
      do 73251 i=1,jmax
      do 73241 j=i+1,jmax
      covinv(i,j)=covinv(j,i)
73241 continue
73251 continue
       if (iswitch_verbose .eq. 1) then

       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' finish inversion', del_sec
       end if

c
c  check for inverse
c
      icheck=0
      if (icheck .eq. 1) then
      print*,' starting check'
      suminv=0.
      do 77110 i=1,jmax
c      write(77,777)i,(covinv(i,j),j=1,jmax)
      do 77111 j=1,jmax
      dsum=0.
        do 77112 k=1,jmax
        dsum=dsum+covar(i,k)*covinv(k,j)
77112    continue
      if (i .eq. j) suminv=suminv+(1.-dsum)**2
      if (i .ne. j) suminv=suminv+dsum**2
77111  continue
77110  continue
      print*,'inverse check',sqrt(suminv/(float(kmax)*float(kmax)))
      end if

      call model_fit(iswt,ModType,ax,d,dr,chi2)
c  need to do transpose!
c      sumres=0.
c      do 73200 i=1,jmax
c      res(i)=0.
c        do 73213 j=1,jmax
c        res(i)=res(i)+dr(i)*covinv(i,j)*dr(j)
c73213   continue
c      print*,i,res(i)
c      sumres=sumres+res(i)
c73200 continue
      if (iswitch_verbose .eq. 1) 
     &  print*,' Z^t * C^{-1} * Z =', chi2
      end if
c End Standard Covariance


c
c  compute maximum likelyhood estimate
c
c  dett is half the true determainant determinant is
c      fmle=(-1.0*dett-0.5*sumres)
c     &      -0.5*float(jmax)*alog(2.0*3.14159265)
      fmle=(-1.0*dett-0.5*chi2)
     &      -0.5*float(jmax-nmiss)*dlog(2.0d+0*3.14159265d+0)
c       print*,"dett,sumres,jmax,nmiss,ixx",dett,sumres,jmax,nmiss,ixx
c      print*,' ln of max likelihood=', fmle
       timex=secnds(0.0)
       del_sec=timex-time0
c      timeZ=secnds(0.0)
c      print*,' MLE done',timeZ-time0

      if (sig4 .ne. sig4_orig) sig4=sig4_orig
      if (sig1 .ne. sig1_orig) sig1=sig1_orig
      if (chi2 .gt. 999999999. ) fmle=-999999999.
c      if (sumres .gt. 999999999. ) sumres=999999999.
      if (fmle .lt. -9999999.) fmle=-9999999.
      if (chi2 .lt. 0.) fmle=-9999999.
      write(6,7491)sig1,sig2,sig3,sig4,sig5,sig6,sig7,dett,
     &   chi2,fmle,del_sec
      write(12,7490)sig1,sig2,sig3,sig4,sig5,sig6,sig7,dett,
     & chi2,fmle
c7490  format(7f13.4,2x,e10.3,2f13.3)
c7491  format(7f13.4,2x,e10.3,3f13.3)
7490  format(7f13.4,2x,3f13.3)
7491  format(7f13.4,2x,4f13.3)
c  print out the time constants if 'float'
      if (nexp .ne. 0) write(6,7493)(bexp(k),k=1,nexp)
7493  format( 10(' tau= ',e10.2, '  '))



       if( iswithc_verbose .eq. 1) 
     &  print*,' finish MLE', del_sec
c      print*,sig4_orig,sig1_orig
      if (sig4 .ne. sig4_orig) sig4=sig4_orig
      if (sig1 .ne. sig1_orig) sig1=sig1_orig
c      print*,sig4_orig,sig1_orig
c
c  rescale noise amplitudes to optimize MLE
c
      nscale=0
      if (nscale .eq. 1) then
      if ((fmle .ne. -9999999.) .or. (chi2 .ne. 999999999.)) then
        scale=sqrt(chi2/float(jmax-nmiss))
        if ((scale .lt. 0.99) .or. (scale .gt. 1.01)) then
c        print*,'scale=',scale
        sig1=sig1*scale
        sig2=sig2*scale
        sig5=sig5*scale
        sig6=sig6*scale
        sumres=chi2/(scale**2)
        dett=dett+float(jmax-nmiss)*alog(scale)
        fmle=(-1.0*dett-0.5*chi2)
     &      -0.5*float(jmax-nmiss)*dlog(2.0d+0*3.14159265d+0)  
        if (chi2 .lt. 0.) fmle=-9999999.
        write(6,7491)sig1,sig2,sig3,sig4,sig5,sig6,sig7,dett,
     &   chi2,fmle,0.0
        write(12,7490)sig1,sig2,sig3,sig4,sig5,sig6,sig7,dett,
     & chi2,fmle   
130     format(90(3x,f15.3,1x,f15.3))
        write(13,130)(x(i),xe(i),i=1,nmod)
        end if
      end if 
      end if
      return
      end

