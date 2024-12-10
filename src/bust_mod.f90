module bust_mod
! subroutines for bust program
!  use iso_fortran_env

  implicit none
  private
 
  public :: i4vec_median, i4vec_frac
  
  contains 
  
      subroutine i4vec_median ( n, a, median )

!*********************************************************************72
!
!! I4VEC_MEDIAN returns the median of an unsorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a ve!tor of I4's.
!
!    Hoare's algorithm is used.  The values of the ve!tor are
!    rearranged by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.

!  Modified:

!    06 July 2009

!  Author:

!    John Burkardt

!  Parameters:

!    Input, integer N, the number of elements of A.

!    Input/output, integer A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.

!    Output, integer MEDIAN, the value of the median of A.

!  MODIFIED BY JL for real MEDIAN and real A
      implicit none

      integer n,i

      real a(n)
      integer k
      real median

      k = ( n + 1 ) / 2
!      print*,n,k,(a(i),i=1,5)

      if (n .ne. 0 ) then
!  BUG FIX BY JOL.... crashes when n=0
      call i4vec_frac ( n, a, k, median )
      else
      median=0
      end if

      return
      end subroutine i4vec_median
      subroutine i4vec_frac ( n, a, k, frac )

!*********************************************************************72

!! I4VEC_FRAC searches for the K-th smallest element in an I4VEC.

!  Discussion:

!    An I4VEC is a vector of I4's.

!    Hoare's algorithm is used.

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    18 July 2010

!  Author:

!    John Burkardt

!  Parameters:

!    Input, integer N, the number of elements of A.

!    Input/output, integer A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.

!    Input, integer K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.

!    Output, integer FRAC, the value of the K-th fractile of A.

      implicit none

      integer n

      real a(n)
      real frac,ix,t
      integer i
      integer iryt
!      integer ix
      integer j
      integer k
      integer left
!      integer t

!  Changed by JOL when n = 0
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

!  Find I so that IX <= A(I).

30        continue

          if ( a(i) .lt. ix ) then
            i = i + 1
            go to 30
          end if

!  Find J so that A(J) <= IX.

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
      end subroutine i4vec_frac 
      subroutine i4vec_heap_a ( n, a )

!*********************************************************************72

! I4VEC_HEAP_A reorders an I4VEC into an ascending heap.

!  Discussion:

!!    An I4VEC is a vector of I4's.

!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).

!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    18 July 2010

!  Author:

!    John Burkardt

!  Reference:

!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.

!  Parameters:

!    Input, integer N, the size of the input array.

!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.

!    MODIFIED by JL to take Real values
      implicit none

      integer n

      real a(n)
      integer i
      integer ifree
      real key
      integer m

!  Only nodes N/2 down to 1 can be "parent" nodes.

      do i = n / 2, 1, -1

!  Copy the value out of the parent node.
!  Position IFREE is now "open".

        key = a(i)
        ifree = i

10      continue

!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)

          m = 2 * ifree

!  Does the first position exist?

          if ( n .lt. m ) then
            go to 20
          end if

!  Does the second position exist?

          if ( m + 1 .le. n ) then

!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.

            if ( a(m+1) .lt. a(m) ) then
              m = m + 1
            end if

          end if

!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.

          if ( key .le. a(m) ) then
            go to 20
          end if

          a(ifree) = a(m)
          ifree = m

        go to 10

!  Once there is no more shifting to do, KEY moves into the free spot.

20      continue

        a(ifree) = key

      end do

      return
      end subroutine i4vec_heap_a 

      subroutine i4vec_sort_heap_d ( n, a )

!*********************************************************************72

!! I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.

!  Discussion:

!!    An I4VEC is a vector of I4's.

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    31 August 2008

!  Author:

!    John Burkardt

!  Reference:

!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.

!  Parameters:

!    Input, integer N, the number of entries in the array.

!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.

!  Modified by JL for real data, not integers

      implicit none

      integer n

      real a(n)
      integer n1

      if ( n .le. 1 ) then
        return
      end if

!  1: Put A into ascending heap form.

      call i4vec_heap_a ( n, a )

!  2: Sort A.

!  The smallest object in the heap is in A(1).
!  Move it to position A(N).

      call i4_swap ( a(1), a(n) )

!  Consider the diminished heap of size N1.

      do n1 = n - 1, 2, -1

!  Restore the heap structure of A(1) through A(N1).

        call i4vec_heap_a ( n1, a )

!  Take the smallest object from A(1) and move it to A(N1).

        call i4_swap ( a(1), a(n1) )

      end do

      return
      end subroutine i4vec_sort_heap_d
      subroutine i4_swap ( i, j )

!*********************************************************************72

!cc I4_SWAP switches two I4's.

!  Licensing:

!    This code is distributed under the GNU LGPL license.

!  Modified:

!    06 January 2006

!  Author:

!    John Burkardt

!  Parameters:

!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.

      implicit none

      real i
      real j
      real k

      k = i
      i = j
      j = k

      return
      end subroutine i4_swap
end module bust_mod