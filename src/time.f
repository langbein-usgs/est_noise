ccc  Subroutines
c   subroutine date2day(yr,mon,date,day)
c   subroutine invcal(day,yr1,mn,date)
c   subroutine jul2day(yr,jul,iday)
c   subroutine inv_jul_time(itime,nyr,jul)
C
C  All calculations are relative to Jan 1, 1960

      subroutine date2day(yr,mon,date,day)
      integer yr,date,day,mon
      dimension mnth(12)
      data (mnth(i),i=1,12)/0,31,59,90,120,151,181,212,243,273,304,334/
c
      yr=yr-1900

      day=0
      day=int(yr-60)*365
      if (yr .le. 60) go to 3
      day= day+1
      if (yr .le. 64) go to 3
      day= day+1
      if (yr .le. 68) go to 3
      day= day+1
      if (yr .le. 72) go to 3
      day= day+1
      if (yr .le. 76) go to 3
      day= day+1
      if (yr .le. 80) go to 3
      day= day+1
      if (yr .le. 84) go to 3
      day=day+1
      if (yr .le. 88) go to 3
      day=day+1
      if (yr .le. 92) go to 3
      day=day+1
      if (yr .le. 96) go to 3
      day=day+1
      if (yr .le. 100) go to 3
      day=day+1     
      if (yr .le. 104) go to 3
      day=day+1     
      if (yr .le. 108) go to 3
      day=day+1     
      if (yr .le. 112) go to 3
      day=day+1     
      if (yr .le. 116) go to 3
      day=day+1 
      if (yr .le. 120) go to 3
      day=day+1     
      if (yr .le. 124) go to 3
      day=day+1     
      if (yr .le. 128) go to 3
      day=day+1     
      if (yr .le. 132) go to 3
      day=day+1         
      if (yr .le. 136) go to 3
      day=day+1
      if (yr .le. 140) go to 3
      day=day+1     
3     day=day+mnth(mon)+date
      if (((yr .eq. 60) .or. (yr .eq. 64) .or. (yr .eq. 68) .or.
     &  (yr .eq. 72) .or. (yr .eq. 76) .or. (yr .eq. 80) .or. 
     & (yr .eq. 84) .or.(yr .eq. 88)
     &  .or. (yr .eq. 92) .or. (yr .eq. 96) .or. (yr .eq. 100)
     & .or. (yr .eq. 104) .or. (yr .eq. 108) .or. (yr .eq. 112)
     & .or. (yr .eq. 116).or. (yr .eq. 120) .or. (yr .eq. 124)
     & .or. (yr .eq. 128).or. (yr .eq. 132) .or. (yr .eq. 136)
     & .or. (yr .eq. 140)
     &  )
     1  .and. (mon .ge. 3)) go to 10
      go to 11
10    day= day+1
      
11    continue
      yr=yr+1900
      return
      end

       subroutine invcal(day,yr1,mn,date)
c
c inverse calender function
c
c  give day number relative to jan 1, 1980
c  output dat yr, month and day
c
       integer day,yr,mn,date,yr1
       integer year(120),month(22),lmnth(22)
c       data (year(i),i=1,20)/0,366,731,1096,1462,1827,2192,2557,2923,
c     1  3288,3653,4018,4384,4749,5114,5479,5845,6210,6575,6940/
       data (month(i),i=1,13)/0,31,59,90,120,151,181,212,243,273,304,
     1                           334,365/
       data (lmnth(i),i=1,13)/0,31,60,91,121,152,182,213,244,274,305,
     1                          335,366/
       icnt=-365
       do 40 j=1,120
         jj=j+2
         xx=float(jj)/4 - (ifix(jj/4.))
         if (xx .ne. 0) add=365
         if (xx .eq. 0) add=366
c         ixx=ifix(jj/4.) - int(jj/4)
c         if (ixx .ne. 0) add=365
c         if (ixx .eq. 0) add=366
         icnt=icnt+add
         year(j)=icnt
40     continue

c select the year
       yr=0
1      continue
       if (day .le. year(yr+1)) go to 5
       yr=yr+1
       go to 1
5      continue
       yr=yr-1
       id=day-year(yr+1)
       yr=yr+60
       yr1=yr+1900
c  leap year sense
       if ((yr .eq. 80) .or. (yr .eq. 84) .or. (yr .eq. 88) .or.
     & (yr .eq. 92) .or. (yr .eq. 96) .or. (yr .eq. 100) .or.
     & (yr .eq. 104) .or. (yr .eq. 108) .or. (yr .eq. 112) .or.
     &   (yr .eq. 116) .or. (yr .eq. 68) .or.
     & (yr .eq. 60) .or. (yr .eq. 64) .or. (yr .eq. 88) .or.
     & (yr .eq. 72) .or. (yr .eq. 76) .or.
     & (yr .eq. 120) .or. (yr .eq. 124) .or. (yr .eq. 128) .or.
     & (yr .eq. 132) .or. (yr .eq. 136) .or. (yr .eq. 140) 
     & ) go to 20
c  select the month
       mn=1
7      continue
       if (id .le. month(mn)) go to 10
       mn=mn+1
       go to 7
10     mn=mn-1
       date=id-month(mn)
       return
c  select month for the leap year
20     continue
       mn=1
27     continue
       if (id .le. lmnth(mn)) go to 30
       mn=mn+1
       go to 27
30     mn=mn-1
       date=id-lmnth(mn)
       return
       end


      subroutine jul2day(yr,jul,iday)
      integer yr, jul, iday
      yr=yr-1900
      iday=0
      if (yr .gt. 60 ) iday=iday+1
      if (yr .gt. 64 ) iday=iday+1
      if (yr .gt. 68 ) iday=iday+1
      if (yr .gt. 72 ) iday=iday+1
      if (yr .gt. 76 ) iday=iday+1
      if (yr .gt. 80 ) iday=iday+1
      if (yr .gt. 84 ) iday=iday+1
      if (yr .gt. 88 ) iday=iday+1
      if (yr .gt. 92 ) iday=iday+1
      if (yr .gt. 96 ) iday=iday+1
      if (yr .gt. 100 ) iday=iday+1
      if (yr .gt. 104 ) iday=iday+1
      if (yr .gt. 108 ) iday=iday+1
      if (yr .gt. 112 ) iday=iday+1
      if (yr .gt. 116 ) iday=iday+1
      if (yr .gt. 120 ) iday=iday+1
      if (yr .gt. 124 ) iday=iday+1
      if (yr .gt. 128 ) iday=iday+1
      if (yr .gt. 132 ) iday=iday+1
      if (yr .gt. 136 ) iday=iday+1
      if (yr .gt. 140 ) iday=iday+1

      iday=iday+(yr-60)*365 + jul
      yr=yr+1900
      return
      end
      subroutine inv_jul_time(itime,nyr,jul)
c converts day number since jan 1, 1960 back to julian day and year
c
c  itime is number of days since Jan 1, 1960
      dimension nyear(120),nleap(120)
      nx=120
        do 1 i=1,nx
        nleap(i)=1+(i-1)*4
1       continue
        j=1
        ncnt=0
        do 2 i=1,nx
        if (i .eq. nleap(j)) then
           ncnt=ncnt+366
           j=j+1
        else
          ncnt=ncnt+365
        end if
        nyear(i)=ncnt
c        print*,i,nyear(i)
2       continue
      k=1
      nyr=1960
c      if (itime .lt. nyear(1)) then
      if (itime .le. nyear(1)) then
        jul=itime
        return
      end if
      k=2
      nyr=nyr+1
3     continue
      if ((itime .gt. nyear(k-1)) .and. (itime .le. nyear(k))) then
        jul=itime-nyear(k-1)
        return
      else
        k=k+1
        nyr=nyr+1
        go to 3
      end if
      end

