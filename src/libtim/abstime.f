c
c
      real*8 function abstime(jyr,jmo,jdy,jhr,jmn,sec)
c
c   Purpose:
c     To convert  a date and time specification  in seconds starting  
c     from 1900 1 1 0 0 0.0.    
c
c   WARNINGS: - since  the routine does not account  for the secular
c               leap year, its validity extends only to 2099.
c             - the routine  does not have a control  on the century
c               and for instance it assumes 1979 to be 79.
c             - if a date is not defined  (ie. it looks like 0 0 0 0
c               0 0.0) the routine gives back 0.0.
c             - if the date is expressed as day of the year, enter 0
c               for the month.
c
c   Arguments:
c     jyr        year.
c     jmo        month. (0 - if day of the year)
c     jdy        day.    
c     jhr        hour.
c     jmn        minute.
c     sec        second.
c     abstime    time in seconds from 1900 1 1 0 0 0.0.
c
c
      real*8 sec
      dimension n(0:12),l(0:12)
      data n/1,0,31,59,90,120,151,181,212,243,273,304,334/
      data l/1,0,31,60,91,121,152,182,213,244,274,305,335/
c
c   Transform date and time in seconds from 1900 1 1 0 0 0.0.
c
      if (jyr.ge.1900) then
         jyear=jyr-1900
      else
         jyear=jyr
      endif
      if (jyear.eq.0) then
         n0=0
         lc=0
         lcdiff=0
      else
         n0=1
         lc=(jyear-1)/4
         lcdiff=(jyear-1)-(lc*4)
      endif
      if (lcdiff.eq.3) then
         jdd=l(jmo)+jdy-1
      else
         jdd=n(jmo)+jdy-1
      endif
      abstime=lc*126230400.0d+00+(lcdiff+n0)*31536000.0d+00
     *        +jdd*86400.0d+00+jhr*3600.0d+00+jmn*60.0d+00+sec
      if (jmo.eq.0.and.jdy.ne.0) abstime=abstime-86400d+00
c
c
      return
      end


