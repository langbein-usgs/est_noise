      subroutine pow_law_cov(corran,jmax,sig1,sig2,sig3,sig4,covar,
     & mmax,md,t_small)
c  Power Law covariance
      double precision covar(md,md),t_small
      double precision dcr(mmax), dsum,corran(md),sumlast
     &  ,dmax,dmin,dmin0,dcmax,dcmin,dsmall,dtemp
      integer irow(md)

      do 1 i=1,jmax
        do 2 j=1,i
        covar(i,j)=0.
        covar(j,i)=0.
2       continue
1     continue




c   power law covariance
c assumes corran(1) .eq. 0
c  colored noise for exponent
c
c USE FRACTIONAL DIFFERENCE
c       print*,'  jmax=', jmax,' corran(jmax)= ',corran(jmax)-corran(1)
       kmax=int((corran(jmax)-corran(1))/t_small + 0.5)+1
       if (kmax .gt. mmax) then
         print*,' kmax = ',kmax,' exceeds mmax of ',mmax
         stop
       end if
       time0=secnds(0.0)
c  modify to fix GGM problem (Apr 2021)
       kkmax=kmax
       irowOffset=0
       if (sig4 .ne. 0.0) then

         kkmax=kmax+ifix(3*1/(sngl(t_small)*sig4))
         irowOffset=kkmax-kmax
c         print*,' kkmax=', kkmax, 'irowOffset=',irowOffset
       end if

       call  frac_diff(dcr,sig3,sig4,sngl(t_small),kkmax,mmax)


c  Bug fix;  noted that standard error of nominal value was incorrect for
c    kswitch0=0; okay for kswitch0=1;  Problem rooted in that first row (or column)
c    of covariance matrix was not all zeros.

      kswitch0=1
c      kswitch0=0
       if ( kswitch0 .ne. 0 ) then
       do 881 i=mmax,2,-1
       dcr(i)=dcr(i-1)
881    continue
       dcr(1)=0
       end if

  
       const=(sig2**2)*(t_small**(sig3/2.0))


c  Original


      if ( kswitch0 .eq. 0) then

C original

       do 7006 i=1,jmax
       imax=int((corran(i)-corran(1))/t_small+0.5)+1

       sumlast=0.
       do 7005 k=1,i
         dsum=0.0
         kmax=int((corran(k)-corran(1))/t_small+0.5)+1


         do 7004  kk=1,kmax


         dsum=dsum+dcr(imax+1-kk)*dcr(kmax+1-kk)
7004     continue


         covar(i,k)=sngl(dsum)*const
         if (i .eq. k) covar(i,k)=covar(i,k)+sig1**2
7005   continue

7006   continue

       else 
c  Modified
       imak=1
       if ( imak .eq. 0) then
         do 70061 i=2,jmax
         imax=int((corran(i)-corran(1))/t_small+0.5)+1

c       if (i .eq. jmax) print*,'jmax and imax= ',jmax,imax
         sumlast=0.
         do 70051 k=2,i
           dsum=0.0
           kmax=int((corran(k)-corran(1))/t_small+0.5)+1
           kkmax=kmax
c         print*,k,kmax
           do 70041  kk=1,kmax


           dsum=dsum+dcr(imax+1-kk)*dcr(kmax+1-kk)
70041       continue

c


           covar(i,k)=sngl(dsum)*const
           if (i .eq. k) covar(i,k)=covar(i,k)+sig1**2
70051     continue

70061     continue
        covar(1,1)=sig1**2
c  modification May 2006
c        covar(1,1)=covar(1,1)+const
        end if
        if (imak .eq. 1)  then
c          do i=1,irowOffset
c             irow(i)=0
c          end do
          do 2011i=1,jmax
            ix=int((corran(i)-corran(1))/t_small+0.5)+1

            irow(i)=ix+irowOffset
2011      continue
         imax=int((corran(jmax)-corran(1))/t_small+0.5)+1+irowOffset


c         print*,' calling makeMatrix'
c         print*,'imax,jamx,mmax,md', imax,jmax,mmax,md
c         print*,'corran', corran(jmax),corran(1),t_small
c         do 3891 i=1,jmax
c           write(56,*)i,irow(i),dcr(irow(i))
c3891     continue
c         write(56,*)' end'
c         do 3892 i=1,imax
c           write(57,*)i,dcr(i)
c3892     continue
c         write(57,*)' end'

         call makeMatrix(irow,imax,(dcr),covar,jmax,0,mmax,md)
            

         do 2801 i=1,jmax
           do 2802 j=1,i
           covar(i,j)=covar(i,j)*const
           if (i .eq. j) covar(i,j)=covar(i,j)+sig1**2
           covar(j,i)=covar(i,j)
c           ccovar(i,j)=covar(i,j)
c           ccovar(j,i)=covar(i,j)
2802       continue
2801     continue

        end if


      end if


      return
      end
      subroutine frac_diff(crr,expon,alpha,ts,npts,nmax)
c use Hosking routine to create fractional differences
c  crr are the fractional difference
c  expon is the exponent
c  alpha is the GM freq in radians/yr
c  ts is the sampling interval in years
c  npts is length of time series
c  nmax is dimension of crr
      double precision crr(nmax),small
      small=1.0d-10
      ddexp=0.5*expon
         do 39 j=1,npts+100
         if (j .eq. 1) fracti=1.0
         if (j .eq. 2) then
            fracti=ddexp
            fracti_last=fracti
         end if
         if (j .ge. 3) then
            fracti=fracti_last*(ddexp+float(j-2))/float(j-1)
            fracti_last=fracti
         end if

         crr(j)=fracti*dexp(-dble(alpha)*ts*float(j-1))
         if (crr(j) .lt. small) crr(j)=small

     
39       continue
      return
      end
      subroutine band_pass_cov(corran,jmax,t_small,flow,fhigh,npole,
     &   mmax,md,sig5,covar)

       double precision corran(md),covar(md,md),cr(mmax),t_small
       len=int((corran(jmax)-corran(1))*365.25+0.5)
       print*,' BP filt called',jmax,t_small,flow,fhigh,
     &  npole,corran(jmax),sig5,mmax,md
      time0=secnds(0.0)
      do 1 i=1,jmax
        do 2 j=1,i
        covar(i,j)=0.
        covar(j,i)=0.
2       continue
1     continue

       call band_pass_filt(t_small,flow,fhigh,npole,mmax,len,cr,sig5)
       print*,' BP filt 2',jmax,t_small,flow,fhigh,
     &  npole,corran(jmax),sig5,mmax,md

       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' form Band Pass filter', del_sec


       do 9056 i=1,jmax
       imax=int((corran(i)-corran(1))/t_small+0.5)+1
       sumlast=0.
       do 9055 k=1,i
         dsum=0.0
         kmax=int((corran(k)-corran(1))/t_small+0.5)+1
         do 9054  kk=1,kmax

         dsum=dsum+cr(imax+1-kk)*cr(kmax+1-kk)
9054     continue
         sumlast=dsum
         covar(i,k)=dsum
9055   continue
9056   continue
c       print*,' Covariance for BP-filtered noise constructed'

       timex=secnds(0.0)
       del_sec=timex-time0
       print*,' form Band Pass filter Covariance', del_sec

      return
      end

      subroutine band_pass_filt(ts,fl,fh,np,nmax,npts,Fi,amp)
c create a bandpass filter in the time domain
c  define filter in frequency domain then take inverse fft of filter
c    tranfer function to get time domain impulse responce
c  ts is sampling interval in years
c  fl and fh are low and high freq break pts in c/yr
c  np is number of poles 1 to 4
c  nmax is the dimension
c  npts is the number of points in time series
c  F on output is the impulse responce
c  amp is the amplitude of filter.
c
c  Modified May 19, 2015 to use fftpack from netlib.org instead of
c   realft from Numerical Recipes

c      dimension H(nmax)
c      double complex hc(nmax),h1,h2,h3,h4
c      double precision Fi(nmax),ts, wsave(262159),a(nmax),b(nmax)
c      dimension H(32768)
c      double complex hc(32768),h1,h2,h3,h4
c      double precision Fi(nmax),ts,wsave(196623),a(32768),b(32768)
cc  swave dim 3*n+15  196623
      double complex hc(32768),h1,h2,h3,h4
      double precision Fi(nmax),ts,wsave(98319),a(32768),b(32768)
     &, sumF
c  swave dim 3*n+15  196623


      fh=fh*ts
      fl=fl*ts
      fn=0.5
      print*,'In BP_filt v2',ts,fl,fh,np,nmax,npts,amp
c  fix bug;
      len=npts+100
      nf=int(len/2)
      nf=int((log(float(len/2))/log(2.))+0)
      nf=nf+1
      nf=2**nf
      print*,' Number of frequencies ',nf
      len=2*nf
      print*,' length of time series for band-pass noise is ',len
c      print*,ts,fl,fh,np,amp
      hc(1)=dcmplx(0.0d+0,0.0d+0)

      peak=-99999.
      do 50 i=1,nf
        f=(1/float(len))*(i)
        h1=dcmplx(1.0,0.0)
        h2=dcmplx(0.,f/fl)
        h3=dcmplx(1.0,f/fl)
        h4=dcmplx(1.,f/fh)
        hc(i)=h1*h2/(h3*h4)
        if (np .eq. 2) hc(i)=hc(i)*hc(i)
        if (np .eq. 3) hc(i)=hc(i)*hc(i)*hc(i)
        if (np .eq. 4) hc(i)=hc(i)*hc(i)*hc(i)*hc(i)
        if (np .eq. 5) hc(i)=hc(i)*hc(i)*hc(i)*hc(i)*hc(i)
        if (np .eq. 6) hc(i)=hc(i)*hc(i)*hc(i)*hc(i)*hc(i)*hc(i)
        a(i)=real(hc(i))
        b(i)=-1.0*aimag(hc(i))

        if (sqrt(a(i)**2+b(i)**2) .gt. peak) then
           peak=sqrt(a(i)**2+b(i)**2)
           f_peak=f/ts
        end if
50    continue
c  Initialize
      call dzffti(len,wsave)


      call dzfftb(len,Fi,0.0d+0,a,b,wsave)
      do 60 i=1,len
      Fi(i)=((amp/peak)*Fi(i)*2.0/float(len))*sqrt(0.5/ts)
60    continue
      fh=fh/ts
      fl=fl/ts
c
c  smooth Fi output using npoint smoother
c    gets rid of tiny oscilations in invFFT of bp filter 
c
c      do i=1,len
c        wsave(i)=Fi(i)
c      end do
c      Fi(2)=(wsave(1)+wsave(2)+wsave(3))/3.0
c      Fi(len-1)=(wsave(len-2)+wsave(len-1)+wsave(len))/3.0
c      do i=3,len-2
c        Fi(i)=(wsave(i-2)+wsave(i-1)+wsave(i)+wsave(i+1)+wsave(i+2))/5.0
c      end do
c
c  smooth Fi output using npoint smoother
c

      npoint=3
      nptt=int(npoint/2)
      do i=1,len
c        if (abs(Fi(i)) .lt. 1.0e-6) Fi(i)=0.0
        wsave(i)=Fi(i)
      end do
      do i=1,len
        sumF=0.0
        nsum=0
c        if ((i .gt. nptt) .and. (i .lt. len-nptt)) then
        do j=1,npoint
          k=i-nptt+j-1
          if (( k .ge. 1) .and. (k .le. len)) then
            sumF=sumF+wsave(k)
            nsum=nsum+1
          end if
        end do
        Fi(i)=sumF/float(nsum)
c        write(46,460)i,wsave(i),Fi(i),wsave(i)-Fi(i),nsum
c460     format(i8,1x,3d15.6,i6)
c        end if
        end do
c        do i=1,len
c          Fi(i)=wsave(i)
c        end do
      return
      end
      subroutine makeMatrix(irow,max,f,covar,nobs,iflip,mmax,nmax)
c  make covariance matrix, covar, from the filter function,f
c  Initially, builds covariance matrix on the assumption unformally sampled observations
c    with no missing data
c  Then, selects rows and columns specified by irow to build actual covariance matrix
c
c     irow  vector contain row numbers (or column numbers) to select out of complete covariance matrxi
c     max  total size of covariance matrix (and f) assuming no missing data
c     f  filter function (length max)
c     covar  covariance matrix
c     nobs  number of data
c     nmax  dimension of covar
c     iflip  if 0; for calculate covariance; if 1--calculating inverse covariance
c     iflip 2  is the CORRECT way of computing the so-called covariance matrix given finv
c
c      Normally, if the H matrix made up of filter coefficents of f (See Langbein 2004, eq 8)
c       Then the cov=H*H'    <------------iflip=0
c     But, if we're given the inverse filter coefficents, finv
c       Then the cov^-1 = H'*H  <---------- iflip=2
c
c      dimension covar(12687,12687), scr(12687,12687),f(12687)
c      dimension covar(12687,12687), scr(12687,12687),f(mmax)
c      double precision  scr(max,max),f(mmax)
c  max dim of scr seems to be 18001,18001 -- any larger, the code crashes immediately

c      double precision covar(12687,12687), scr(max,max),f(mmax)
c      double precision covar(nmax,nmax), scr(max,max),f(mmax)
c      double precision covar(nmax,nmax), scr(32768,32768),f(mmax)
c      double precision covar(12687,12687), scr(32768,32768),f(32768)
c  compiles and run with  intel/mac but doesn't compile with gfortran/mac
c      double precision covar(nmax,nmax), scr(32768,32768),f(mmax)
c compiles and runs with gfortan/mac  compiles but doesn't run with intel/mac
c      double precision covar(nmax,nmax), scr(mmax,mmax),f(mmax)

c  compiles and run with  intel/mac but doesn't compile with gfortran/mac
      double precision covar(nmax,nmax), scr(32758,32768),f(mmax)

c      integer irow(12687)
      integer irow(nmax)
      double precision sum,tmp

c      print*,' max nobs, mmax nmax', max, nobs, mmax, nmax
      if (iflip .lt. 2 ) then
      if (nobs .ne. max ) then
c  do this when there are missing data
      do 10 i=1,max
        imax=i
        do 15 k=1,i
        sum=0.0
        kmax=k
c         if (( i .ge. 1) .or. (k .ge. 1 )) then
        if (( i .eq. 1) .or. (k .eq. 1 )) then
          sum1=0.0
          sum2=0.0
          kmin=kmax
        else
          sum=scr(i-1,k-1)
          kmin=1
        end if
        do 203 kk=1,kmin

          sum=sum+f(imax+1-kk)*f(kmax+1-kk)
203     continue
        scr(i,k)=sum
        scr(k,i)=sum

15      continue

10    continue
c
c  select row/columns
c
      do 20 i=1,nobs
        krow=irow(i)
        do 25 j=1,i
        kcol=irow(j)
        covar(i,j)=scr(krow,kcol)
        covar(j,i)=scr(kcol,krow)
25      continue
20    continue


      else
c  when there are no missing data, do this....  

      do 11 i=1,max
        imax=i
        do 16 k=1,i
        sum=0.0
        kmax=k
c         if (( i .ge. 1) .or. (k .ge. 1 )) then
        if (( i .eq. 1) .or. (k .eq. 1 )) then
          sum1=0.0
          sum2=0.0
          kmin=kmax
        else
          sum=covar(i-1,k-1)
          kmin=1
        end if
        do 204 kk=1,kmin

          sum=sum+f(imax+1-kk)*f(kmax+1-kk)
204     continue
        covar(i,k)=sum
        covar(k,i)=sum

16      continue

11    continue
      end if
c      do 3331 i=1,100
c        scr(i,1)=0.d+0
c3331  continue
c      do 3332 i=1,100
c        write(99,*)f(i)
c        write(98,980)(covar(i,j),j=1,1000)
c        write(97,980)(f(k),k=i,1,-1),(scr(k,1),k=i+1,100)
c980    format(1000e20.12)
c3332  continue
c      if (nobs .eq. 1000 ) stop
c      print*,"In makeMatrix f=",(i,f(i),i=1,10)
c      do 3333 i=1,10
c3333  print*,(i,j,covar(i,j),j=1,10)   
c
c   flip matrix for inverse covariance
c
c  NOTE, I'm not sure that this is required; tests suggest with iflip=0,
c   that I get an inverse
      if (iflip .eq. 1 ) then
        nobs2=int(nobs/2)
        idif=nobs-2*nobs2
c   idif=0 then, nobs is even; idif=1 then nobs is odd
        do 30 i=1,nobs2
          do 31 j=1,nobs
            print*,i,j,nobs+1-i,nobs+1-j
          tmp=covar(nobs+1-i,nobs+1-j)
          covar(nobs+1-i,nobs+1-j)=covar(j,i)
          covar(j,i)=tmp
31        continue
30      continue
        if (idif .ne. 0) then
c          i=nobs2+1
          do 33 i=1,nobs2+1
          do 32 j=1,nobs
            print*,i,j,nobs+1-i,nobs+1-j
            tmp=covar(nobs+1-i,nobs+1-j)
            covar(nobs+1-i,nobs+1-j)=covar(j,i)
            covar(j,i)=tmp
32        continue 
33        continue           
        end if

      end if
c  end iflip .ne. 2
      end if
c
c   Do this for inverse covariance matrix given the inverse filter
c
c    NOTE that nobs must equal max
      if (iflip .eq. 2 ) then
      do 41 i=1,max
c      do 41 i=max,1,-1
        imax=i
        do 42 k=1,i
c        do 42 k=i,1,-1
        sum=0.0
        kmax=k
c         if (( i .ge. 1) .or. (k .ge. 1 )) then
        if (( i .eq. 1) .or. (k .eq. 1 )) then
          kmin=kmax
        else
          sum=covar(i-1,k-1)
          kmin=1
        end if
        do 43 kk=1,kmin
c        do 43 kk=kmin,1,-1
c         print*,i,k,kk,imax+1-kk,kmax+1-kk
          sum=sum+f(imax+1-kk)*f(kmax+1-kk)
43      continue
        covar(i,k)=sum
        covar(k,i)=sum

42      continue

41    continue
c
c  flip around the axis from lower left to upper right....
      do 50 i=1,max
        do 51 j=1,max
        scr(max+1-j,max+1-i)=covar(i,j)
51      continue
50    continue
      do 52 i=1,max
      do 53 j=1,i
      covar(i,j)=scr(i,j)
      covar(j,i)=scr(i,j)
53    continue
52    continue
c50    continue
      end if
      return
      end
      subroutine convolv1(npts,H,din,dout,maxtime,maxdata)
c
c  H is filter
c  din and dout are input and output time series
      double precision din(maxdata),dout(maxdata)
      double precision H(maxtime)
c      print*,npts

      do 52 i=1,npts
        dout(i)=0.0
        do 53 j=1,i
        dout(i)=dout(i)+din(j)*H(i-j+1)
53      continue
52    continue
c      do i=1,npts
c        write(6,100)i,H(i),din(i),dout(i)
c      end do
c100   format(i4,f10.5,1x,f10.5,f10.5)
      return
      end

