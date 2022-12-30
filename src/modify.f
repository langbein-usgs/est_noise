      subroutine modify_A(A,t)

c  modifies A matrix to include exponential terms
      double precision A(max_data,max_mod)
      dimension x(82),e(82),bexp(82)
      double precision texp(82),t_start,t(max_data),t_stop
      character*7 exp_choice(10),exp_type(10)
c   double precision

      common /ModFit1a/t_start,t_stop,texp
c  single precision
      common /ModFit2/x,e,bexp
      common /ModFit2a/expmax
c integer
      common /ModFit3/max_data,max_mod,ic,nmod,n_exp,nmiss,irowmiss
c  character
      common /ModFit4/ exp_choice,exp_type
c
c  for each exp_choice equal to "float", add an additional column to A
c
      nexp=n_exp
      tmax=t(1)
      tmin=t(1)
      dtmin=99999.
      ratmax=5.0
      do 20 i=2,ic
        if (t(i) .ge. tmax) tmax=t(i)
        if (t(i) .lt. tmin) tmin=t(i)
        if (t(i)-t(i-1) .lt. dtmin) dtmin=t(i)-t(i-1)
20    continue

      nexp=n_exp

      do 1 kexp=1,nexp
        
        if (exp_choice(kexp) .eq. "float") then
        nmod=nmod+1
        write(89,*)kexp,exp_choice(kexp),texp(kexp),bexp(kexp)
c  test for bexp being too big or too small
        iflagsmall=0
        if (bexp(kexp) .lt.  dtmin/(365.25*ratmax) ) then
           bexp(kexp)= dtmin/(365.25*ratmax) 
           iflagsmall=0
        end if
        if (bexp(kexp) .gt. (tmax-tmin)*ratmax/365.25 ) 
     &        bexp(kexp)= (tmax-tmin)*ratmax/365.25
        if ((bexp(kexp) .gt. (tmax-tmin)*1.0/365.25) .and.
     &         (exp_type(kexp) .eq. "m") ) 
     &        bexp(kexp)= (tmax-tmin)*1.0/365.25
        do 2 j=1,ic
        A(j,nmod)=0.0
        if (exp_type(kexp) .eq. "e") then
        if (t(j) .ge. texp(kexp)-t_start) 
     &    A(j,nmod)=1.0-exp(-(t(j)-(texp(kexp)-t_start))
     &           /(365.25*bexp(kexp)))
         end if
        if (exp_type(kexp) .eq. "o") then
        if ((t(j) .gt. texp(kexp)-t_start) .and. (iflagsmall .eq. 0)) 
     &    A(j,nmod)=alog10(abs(bexp(kexp)) +
     &     (sngl(t(j)-(texp(kexp)-t_start))/365.25))
        if ((t(j) .gt. texp(kexp)-t_start) .and. iflagsmall .eq. 1) 
     &    A(j,nmod)=alog10((sngl(t(j)-(texp(kexp)-t_start))/365.25))
c     &             -alog10(abs(bexp(k)))
         end if
        if (exp_type(kexp) .eq. "m") then
        if ((t(j) .gt. texp(kexp)-t_start) .and. (iflagsmall .eq. 0))
     &    A(j,nmod)=alog10( 1.0 +
     &     (sngl(t(j)-(texp(kexp)-t_start))/365.25)/abs(bexp(kexp)) )
        if ((t(j) .gt. texp(kexp)-t_start) .and. (iflagsmall .eq. 1))
     &    A(j,nmod)=alog10((sngl(t(j)-(texp(kexp)-t_start))/365.25))
     &             -alog10(abs(bexp(k)))
         end if
 
c        write(49,490)j,t(j),(texp(kexp)-t_start),A(j,nmod)
c490     format(i8,3f15.5)
2       continue

        end if

1     continue
      return
      end

