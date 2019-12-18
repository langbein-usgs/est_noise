c attempts to remove observations that are either too large or too small
c
c  specify input time series
c
c  compute a running median (not running mean)
c    after specifying the period
c  
c  examine the residuals about the median
c
c   identify observations greater or less than a threshold around the median
c   eliminate those observations
c
c  difference between bust_2 and bust_3 is that new verision allows for 2 different
c     date formats
c 
c  Bust_4 uses  more efficient algorithms to take medians and/or sort data
c
c  Modified bust_4 (now _5) to:
c  1) remove sorting routines from Press et al Numerical recipes; substituting
c      routines written by John Burkardt; https://people.sc.fsu.edu/~jburkardt/f77_src/i4lib/i4lib.f
c  2) Provide more variety in data formats
c  3) Get time subroutines from time.f
c
      character*20 ifile,ofile
      character*80 string
      character*3 net
      dimension d(9000000),e(9000000),dsrt(9000000),r(9000000)
      integer year(9000000)
      double precision fjul,fday,tsam_day,time0,time,day(9000000),
     &  dec_timed
      fmiss=9999999.
      
      nmax=9000000
      do 1 i=1,nmax
      r(i)=fmiss
1     d(i)=fmiss
      print*,' Input the data format'
      print*,' otr=data'
      print*,' otr   format of year day_of_year obs err'

      print*,' otd   format of YearMnDa  obs err'

      print*,' otx   format of year mo da obs err'

      print*,' mjd   Modified Julian day, obs, err'
      print*,'For all of the above, the da and doy are ',
     &  'double precision'
      print*,'    which allows decimal days to accomodate seconds'
      print*,' '
      print*,' gmt  GMT format for time'
      print*,'        year-mo-daThr:mn:secs.x  obs  err'

      read(5,*)net
105   format(a3)
      if (net .eq. '1') net='otr'
      print*,'net is ', net
      print*,' Input file to examine'
      print*,'   input format is year, julian day, data, data_error'
      read(5,100)ifile
100   format(a20)
      open(1,file=ifile)
      print*,' Input name of output file'
      read(5,100)ofile
      open(4,file=ofile)
      open(3,file="reject.out")
      print*,' Input sampling interval in days'
      print*,'   1--hour = ',1./24.d+0
      print*,'   10--minutes =',10./(24.d+0 * 60.0d+0)
      print*,'   1--minute = ',1.0/(24.d+0 * 60.0d+0)
      print*,'   1--second =', 1.0/(24.d+0 * 3600.0d+0)
      read(5,*)tsam_day
      n=0
      ic=0
2     continue
        read(1,106,end=3) string
106     format(a80)
        call GetData(net,string,time,obs,err)

        if (n .eq. 0) time0=time
c        if (n .eq. 0) then
c           print*,yr,fjul,time
c        end if
        n=int(((time-time0)/tsam_day)+0.5)+1
        ic=ic+1
        day(n)=time
        if (n .gt. nmax) then
          print*,' Number of observation exceeds dimensions'
          print*,' n=',n,' nmax= ',nmax
          print*,' time0=',time0
          print*,' last record read is',yr,fjul,obs
          print*,'  time=', time
          stop
        end if

c        print*,n,n+1,year(n+1),fjul,day(n+1)




c      write(2,*)n,time(n)

      d(n)=obs
      dsrt(ic)=obs
      e(n)=err
c        write(10,*)n,yr,fjul,obs,time,year(n),day(n),d(n)
c      print*,year,fjul,time,ic,n,time-time0,(time-time0)/tsam_day

      if (n .ge. nmax) then
         print*,' Number of observations exceed dimensions'
         print*,'  nmax=',nmax,'  n=',n
         print*,' time0=', time0
         print*,n,yr,fjul,obs,time,year(n),day(n),d(n)

         stop
      end if
      
      go to 2
3     continue
     
      print*,' number of data',ic
      print*,' Total number of days ',n*tsam_day



c  remove median and a secular rate
      if (ic .gt. 0 ) then
        call  i4vec_median (ic, dsrt, xmed )
      else
        xmed=0.0
      end if
      print*,' the median is: ', xmed
      do 12 i=1,n
      if (d(i) .ne. fmiss) r(i)=d(i)-xmed
12    continue
      tstart=0

      print*,' Input the window length in days to compute running ',
     & 'median'
      read(5,*)per
      nwind=int(per/tsam_day)
      print*,'  Number of points in window', nwind
      n_tot_wind=int(n/nwind)
      print*,'  Total number of windows is', n_tot_wind+1
c
c  compute medians of windows of data with length perio 
c
c      istart=1
c      istop=istart+nwind-1
      do 20 k=1,n_tot_wind
      istart=k*nwind-nwind+1
      istop=k*nwind
      ix=0
        do 21 i=istart,istop

           if (r(i) .ne. fmiss) then
           ix=ix+1
           dsrt(ix)=r(i)
           end if
21      continue
c        xmed=fmedian(ix,dsrt)
c        print*,istart,istop,ix
        if (ix .gt. 0) then
          call  i4vec_median (ix, dsrt, xmed )
        else
          xmed=0.0
        end if
        print*,istart,istop,ix,xmed
c   remove running median
        do 22 i=istart,istop
        if (r(i) .ne. fmiss) r(i)=r(i)-xmed
22      continue
20    continue
c  same as above but for the last window
        istart=istop+1
        istop=n
        ix=0
        do 23 i=istart,istop
          if (r(i) .ne. fmiss) then
          ix=ix+1
          dsrt(ix)=r(i)
          end if
23      continue
c        xmed=fmedian(ix,dsrt)
        if (ix .gt. 0) then
          call  i4vec_median (ix, dsrt, xmed )
        else
          xmed=0.0
        end if
        print*,istart,istop,ix,xmed
        do 24 i=istart,istop
        if (r(i) .ne. fmiss) r(i)=r(i)-xmed
24      continue
c
c   sort the data to get statistics on data
c
      sum=0.0
      ix=0
      do 30 i=1,n
      if (r(i) .ne. fmiss) then
        ix=ix+1
        dsrt(ix)=r(i)
        sum=sum+r(i)**2
c        write(98,*)ix,dsrt(ix)
      end if
30    continue
      print*,"  Number of missing observations is ",n-ix
c      call i4vec_heap_a ( ix, dsrt )
c      if (ix .gt. 0 ) then
c        call i4vec_sort_heap_d ( ix, dsrt )
c      end if
c      do 301 i=1,ix
c       write(99,*)i,dsrt(i)
c301   continue
c      call hpsort(ix,dsrt)

c      print*," Minimum value ",1,dsrt(1)
c      print*," Maximum value ",ix,dsrt(ix)
      print*, 'Standard deviation of filtered data is ',
     &   sqrt(sum/float(ix))
      call i4vec_frac ( ix, dsrt, 1, xx1 )
      call i4vec_frac ( ix, dsrt, ix, xx2 )
      print*," Extreme values ",xx1,xx2,
     & " distance = ",xx2-xx1
      n75=int(0.125*ix)
      n25=ix-n75+1
      call i4vec_frac ( ix, dsrt, n75, xx1)
      call i4vec_frac ( ix, dsrt, n25, xx2)
      print*," 75% interval",xx1, xx2,
     & " distance= ",xx2-xx1
      thres0=xx2-xx1
      n90=int(0.05*ix)
      n10=ix-n90+1
      call i4vec_frac ( ix, dsrt, n90, xx1)
      call i4vec_frac ( ix, dsrt, n10, xx2)
      print*," 90% interval",xx1,xx2,
     & " distance= ",xx2-xx1
      n95=int(0.025*ix)      
      n05=ix-n95+1
      if (n95 .lt. 1) n95=1
      if (n05 .gt. ix) n05=ix
      call i4vec_frac ( ix, dsrt, n95, xx1)
      call i4vec_frac ( ix, dsrt, n05, xx2)
      print*," 95% interval",xx1,xx2,
     & " distance= ",xx2-xx1
      n99=int(0.005*ix)
      n01=ix-n99+1
      if (n99 .lt. 1) n99=1
      if (n01 .gt. ix) n01=ix
      call i4vec_frac ( ix, dsrt, n99, xx1)
      call i4vec_frac ( ix, dsrt, n01, xx2)
      print*," 99% interval",xx1,xx2,
     & " distance= ",xx2-xx1
      n999=int(0.0005*ix)
      n001=ix-n999+1
      if (n999 .lt. 1) n999=1
      if (n001 .gt. ix) n001=ix
      call i4vec_frac ( ix, dsrt, n999, xx1)
      call i4vec_frac ( ix, dsrt, n001, xx2)
      print*," 99.9% interval",xx1,xx2,
     & " distance= ",xx2-xx1
      print*, " Input the threshold used to reject data"
      print*,'    use units of IQR'
      read(5,*)thres
      thres=thres*thres0/2.0
      print*,'  The threshold is ',thres
      nbust=0

       do 40 i=1,n
        if (d(i) .ne. fmiss) then
          time=day(i)
          itime=int(time)
          dec_timed=(time)-int(time)
          if (abs(r(i)) .lt. thres) then
           if (net .eq. 'otr' ) then
             call inv_jul_time(itime,nyr,jul)
             write(4,6501)nyr,float(jul)+dec_timed,d(i),e(i)
           end if
           if (net .eq. 'otx' ) then
             call invcal(itime,nyr,mn,idate)
             write(4,6502)nyr,mn,idate+dec_timed,d(i),e(i)
           end if
           if (net .eq. 'otd' ) then
             call invcal(itime,nyr,mn,idate)
             write(4,6503)nyr,mn,idate,dec_timed,d(i),e(i)
           end if 
           if (net .eq. 'gmt' ) then
             call invcal(itime,nyr,mn,idate)
             ihr1=int(24.0*dec_timed)
             imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
             sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
             sec1=sec1*3600.0*24.0
             isec1=int(sec1)
             write(4,6504)nyr,mn,idate,ihr1,imn1,isec1,sec1,d(i),e(i)
           end if
           if (net .eq. 'mjd' ) then
              write(4,6505)dble(itime)+dec_timed+36933.d+0,d(i),e(i)
           end if
          else 
c  output rejected data
          if (net .eq. 'otr' ) then
             call inv_jul_time(itime,nyr,jul)
             write(3,65013)nyr,float(jul)+dec_timed,r(i),d(i),e(i)
           end if
           if (net .eq. 'otx' ) then
             call invcal(itime,nyr,mn,idate)
             write(3,65023)nyr,mn,idate+dec_timed,r(i),d(i),e(i)
           end if
           if (net .eq. 'otd' ) then
             call invcal(itime,nyr,mn,idate)
             write(3,65033)nyr,mn,idate,dec_timed,r(i),d(i),e(i)
           end if 
           if (net .eq. 'gmt' ) then
             call invcal(itime,nyr,mn,idate)
             ihr1=int(24.0*dec_timed)
             imn1=int(24.0*60.0*(dec_timed-dble(ihr1)/24.0))
             sec1=dec_timed-dble(ihr1)/24.0-dble(imn1)/(24.0*60.0)
             sec1=sec1*3600.0*24.0
             isec1=int(sec1)
             write(3,65043)nyr,mn,idate,ihr1,imn1,isec1,sec1,r(i)
     &       ,d(i),e(i)
           end if
           if (net .eq. 'mjd' ) then
              write(3,65053)dble(itime)+dec_timed+36933.d+0,r(i)
     &       ,d(i),e(i)
           end if
c          write(3,*)i,day(i),r(i),d(i),e(i)
          nbust=nbust+1
        end if
        end if      
40     continue
6501  format(1x,i5,1x, f15.9, 1x, f15.2,1x,f8.2)
6502  format(1x,2i5,f15.9, 1x,f15.2, 1x,f8.2)
6503  format(1x,i4,i2.2,i2.2,f10.9, 1x,f15.2,1x,f8.2)
6504  format(1x,i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1, 1x,f15.2,1x,f8.2)
6505  format(1x,f18.9, 1x,f15.2,1x,f8.2)
65013  format(1x,i5,1x, f15.9, 1x, f15.2,1x,f8.2,1x,f8.2)
65023  format(1x,2i5,f15.9, 1x,f15.2, 1x,f8.2,1x,f8.2)
65033  format(1x,i4,i2.2,i2.2,f10.9, 1x,f15.2,1x,f8.2,1x,f8.2)
65043  format(1x,i4,'-',i2.2,'-',i2.2,'T',
     & i2.2,':',i2.2,':',i2.2,f2.1, 1x,f15.2,1x,f8.2,1x,f8.2)
65053  format(1x,f18.9, 1x,f15.2,1x,f8.2,1x,f8.2)




        
        print*," Number of data rejected is ", nbust
c        istart=istop+nwind

       print*," Clean data in ", ofile
       print*,"  Outlyers in file reject.out"



      close (1)
      close (2)
      close (3)
      close (4)
      stop
      end
c  Note --- only i4vec_median and i4vec_frac are used
c   the remaining subroutines below (i4...) aren't needed
      subroutine i4vec_median ( n, a, median )

c*********************************************************************72
c
cc I4VEC_MEDIAN returns the median of an unsorted I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Hoare's algorithm is used.  The values of the vector are
c    rearranged by this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, integer A(N), the array to search.  On output,
c    the order of the elements of A has been somewhat changed.
c
c    Output, integer MEDIAN, the value of the median of A.
c
c  MODIFIED BY JL for real MEDIAN and real A
      implicit none

      integer n,i

      real a(n)
      integer k
      real median

      k = ( n + 1 ) / 2
c      print*,n,k,(a(i),i=1,5)

      if (n .ne. 0 ) then
c  BUG FIX BY JOL.... crashes when n=0
      call i4vec_frac ( n, a, k, median )
      else
      median=0
      end if

      return
      end
      subroutine i4vec_frac ( n, a, k, frac )

c*********************************************************************72
c
cc I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Hoare's algorithm is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, integer A(N), array to search.  On output,
c    the elements of A have been somewhat rearranged.
c
c    Input, integer K, the fractile to be sought.  If K = 1, the
c    minimum entry is sought.  If K = N, the maximum is sought.
c    Other values of K search for the entry which is K-th in size.
c    K must be at least 1, and no greater than N.
c
c    Output, integer FRAC, the value of the K-th fractile of A.
c
      implicit none

      integer n

      real a(n)
      real frac,ix,t
      integer i
      integer iryt
c      integer ix
      integer j
      integer k
      integer left
c      integer t

c  Changed by JOL when n = 0
      if (n .eq. 0 ) then
        frac=0
        return
      end if
      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
        stop 1
      end if

      if ( k .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
        stop 1
      end if

      if ( n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
        stop 1
      end if

      left = 1
      iryt = n

10    continue

        if ( iryt .le. left ) then
          frac = a(k)
          go to 60
        end if

        ix = a(k)
        i = left
        j = iryt

20      continue

          if ( j .lt. i ) then

            if ( j .lt. k ) then
              left = i
            end if

            if ( k .lt. i ) then
              iryt = j
            end if

            go to 50

          end if
c
c  Find I so that IX <= A(I).
c
30        continue

          if ( a(i) .lt. ix ) then
            i = i + 1
            go to 30
          end if
c
c  Find J so that A(J) <= IX.
c
40        continue

          if ( ix .lt. a(j) ) then
            j = j - 1
            go to 40
          end if

          if ( i .le. j ) then

            t    = a(i)
            a(i) = a(j)
            a(j) = t

            i = i + 1
            j = j - 1

          end if

        go to 20

50      continue

      go to 10

60    continue

      return
      end
      subroutine i4vec_heap_a ( n, a )

c*********************************************************************72
c
cc I4VEC_HEAP_A reorders an I4VEC into an ascending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    An ascending heap is an array A with the property that, for every index J,
c    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
c    2*J and 2*J+1 are legal).
c
c                  A(1)
c                /      \
c            A(2)         A(3)
c          /     \        /  \
c      A(4)       A(5)  A(6) A(7)
c      /  \       /   \
c    A(8) A(9) A(10) A(11)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the input array.
c
c    Input/output, integer A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
c    MODIFIED by JL to take Real values
      implicit none

      integer n

      real a(n)
      integer i
      integer ifree
      real key
      integer m
c
c  Only nodes N/2 down to 1 can be "parent" nodes.
c
      do i = n / 2, 1, -1
c
c  Copy the value out of the parent node.
c  Position IFREE is now "open".
c
        key = a(i)
        ifree = i

10      continue
c
c  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
c  IFREE.  (One or both may not exist because they exceed N.)
c
          m = 2 * ifree
c
c  Does the first position exist?
c
          if ( n .lt. m ) then
            go to 20
          end if
c
c  Does the second position exist?
c
          if ( m + 1 .le. n ) then
c
c  If both positions exist, take the smaller of the two values,
c  and update M if necessary.
c
            if ( a(m+1) .lt. a(m) ) then
              m = m + 1
            end if

          end if
c
c  If the small descendant is smaller than KEY, move it up,
c  and update IFREE, the location of the free position, and
c  consider the descendants of THIS position.
c
          if ( key .le. a(m) ) then
            go to 20
          end if

          a(ifree) = a(m)
          ifree = m

        go to 10
c
c  Once there is no more shifting to do, KEY moves into the free spot.
c
20      continue

        a(ifree) = key

      end do

      return
      end

      subroutine i4vec_sort_heap_d ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
c  Modified by JL for real data, not integers
c
      implicit none

      integer n

      real a(n)
      integer n1

      if ( n .le. 1 ) then
        return
      end if
c
c  1: Put A into ascending heap form.
c
      call i4vec_heap_a ( n, a )
c
c  2: Sort A.
c
c  The smallest object in the heap is in A(1).
c  Move it to position A(N).
c
      call i4_swap ( a(1), a(n) )
c
c  Consider the diminished heap of size N1.
c
      do n1 = n - 1, 2, -1
c
c  Restore the heap structure of A(1) through A(N1).
c
        call i4vec_heap_a ( n1, a )
c
c  Take the smallest object from A(1) and move it to A(N1).
c
        call i4_swap ( a(1), a(n1) )

      end do

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      real i
      real j
      real k

      k = i
      i = j
      j = k

      return
      end
