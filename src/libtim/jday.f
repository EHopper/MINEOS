c
c
      integer*4 function jday(jyr,jmo,jdy)
c
c   Purpose:
c     To find out the day of the year (julian day).
c     The routine returns zero if some implausible value is entered.
c
c   Arguments:
c     jyr        year.
c     jmo        month.
c     jdy        day.    
c     jday       day of the year (julian day).
c
c
      logical*4 leap
      integer*4 jyr,jmo,jdy
      integer*4 ilp100,ilp400
      integer*4 n(12),l(12)
      integer*4 jdym(12)
c
      data n/0,31,59,90,120,151,181,212,243,273,304,334/
      data l/0,31,60,91,121,152,182,213,244,274,305,335/
      data jdym/31,29,31,30,31,30,31,31,30,31,30,31/
c
c   Check for absurd values end eventually return.
c
      if (jmo.le.0.or.jmo.ge.13) then
         jday=0
         return
      else if (jdy.le.0.or.jdy.gt.jdym(jmo)) then
         jday=0
         return
      endif
c
c   Check if it is a leap year.
c
      leap=.false.
      ilp=mod(jyr,4)
      if (ilp.eq.0) then
         leap=.true.
         ilp100=mod(jyr,100)
         ilp400=mod(jyr,400)
         if (ilp100.eq.0.and.ilp400.ne.0) leap=.false.
      endif
c
c   Find out the day of the year.
c
      if (leap) then
         jday=l(jmo)+jdy
      else
         jday=n(jmo)+jdy
      endif
c
c   Last check.
c
      if (.not.leap.and.(jmo.eq.2.and.jdy.eq.29)) then
         jday=0
      endif
c
c
      return
      end


