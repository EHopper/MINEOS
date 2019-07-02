c
c
      subroutine time(abstime,jyr,jmo,jdy,jhr,jmn,sec)
c
c   Purpose:
c     To transform  a period of time, expressed as seconds from 1900
c     1 1 0 0 0.0, in a date and time.     
c
c   WARNINGS: - the routine works only up to the year 1999 (2004)
c             - the routine  does not have a control  on the century
c               and for istance it assumes 1979 to be 79. (NO!!!)
c
c   02/02/03 -- modified to return 4-digit year -- above two warnings
c               obsolete.  see comments for syr in catbin program
c               (usr/local/src/Cmt).  Currently up to date through 2004.
c
c   Arguments:
c     abstime    time in seconds from 1900 1 1 0 0 0.0.
c     jyr        year.
c     jmo        month.
c     jdy        day.
c     jhr        hour.
c     jmn        minute.
c     sec        second.
c
c
      implicit real*8 (a-h,o-z)
      parameter (iyrm = 110)
      dimension syr(0:iyrm),smon(13),smol(13)
c
c   The  data statement  of syr  is splitted  because of compliation
c   problems.
c
      data (syr(i),i=0,48)/0.0d+00,
     *0.31536000d+08,0.63072000d+08,0.94608000d+08,0.12614400d+09,
     *0.15776640d+09,0.18930240d+09,0.22083840d+09,0.25237440d+09,
     *0.28399680d+09,0.31553280d+09,0.34706880d+09,0.37860480d+09,
     *0.41022720d+09,0.44176320d+09,0.47329920d+09,0.50483520d+09,
     *0.53645760d+09,0.56799360d+09,0.59952960d+09,0.63106560d+09,
     *0.66268800d+09,0.69422400d+09,0.72576000d+09,0.75729600d+09,
     *0.78891840d+09,0.82045440d+09,0.85199040d+09,0.88352640d+09,
     *0.91514880d+09,0.94668480d+09,0.97822080d+09,0.10097568d+10,
     *0.10413792d+10,0.10729152d+10,0.11044512d+10,0.11359872d+10,
     *0.11676096d+10,0.11991456d+10,0.12306816d+10,0.12622176d+10,
     *0.12938400d+10,0.13253760d+10,0.13569120d+10,0.13884480d+10,
     *0.14200704d+10,0.14516064d+10,0.14831424d+10,0.15146784d+10/
      data (syr(i),i=49,iyrm)/
     *0.15463008d+10,0.15778368d+10,0.16093728d+10,0.16409088d+10,
     *0.16725312d+10,0.17040672d+10,0.17356032d+10,0.17671392d+10,
     *0.17987616d+10,0.18302976d+10,0.18618336d+10,0.18933696d+10,
     *0.19249920d+10,0.19565280d+10,0.19880640d+10,0.20196000d+10,
     *0.20512224d+10,0.20827584d+10,0.21142944d+10,0.21458304d+10,
     *0.21774528d+10,0.22089888d+10,0.22405248d+10,0.22720608d+10,
     *0.23036832d+10,0.23352192d+10,0.23667552d+10,0.23982912d+10,
     *0.24299136d+10,0.24614496d+10,0.24929856d+10,0.25245216d+10,
     *0.25561440d+10,0.25876800d+10,0.26192160d+10,0.26507520d+10,
     *0.26823744d+10,0.27139104d+10,0.27454464d+10,0.27769824d+10,
     *0.28086048d+10,0.28401408d+10,0.28716768d+10,0.29032128d+10,
     *0.29348352d+10,0.29663712d+10,0.29979072d+10,0.30294432d+10,
     *0.30610656d+10,0.30926016d+10,0.31241376d+10,0.31556736d+10,
     *0.31872960d+10,0.32188320d+10,0.32503680d+10,0.32819040d+10,
     *0.33135264d+10,0.33450624d+10,0.33765984d+10,0.34081344d+10,
     *0.34397568d+10,0.1d+11/
      data smon/0.0d+00,0.267840d+07,0.509760d+07,0.777600d+07,
     *0.103680d+08,0.130464d+08,0.156384d+08,0.183168d+08,
     *0.209952d+08,0.235872d+08,0.262656d+08,0.288576d+08,0.1d+09/
      data smol/0.0d+00,0.267840d+07,0.518400d+07,0.786240d+07,
     *0.104544d+08,0.131328d+08,0.157248d+08,0.184032d+08,
     *0.210816d+08,0.236736d+08,0.263520d+08,0.289440d+08,0.1d+09/
      data sdy/0.86400d+05/
      data shr/0.3600d+04/
      data smn/0.60d+02/
c
      abstim=abstime
c
c   Outermost loop to find out the year.
c
      do i=0,iyrm
      if (abstim.ge.syr(i).and.abstim.lt.syr(i+1)) then
      jyr=1900+i
      abstim=abstim-syr(i)
      lc=(i-1)/4
      lcdiff=(i-1)-(lc*4)
      if(lcdiff.eq.3) then
c
c   Second outermost loop to find out the month for a leap year.
c
      do m=1,12
      if(abstim.ge.smol(m).and.abstim.lt.smol(m+1)) then
      jmo=m
      abstim=abstim-smol(m)
c
c   Third outermost loop to find out the day for a leap year.
c
      do l=1,31
      sdymo=(l-1)*sdy
      sdymop1=l*sdy
      if(abstim.ge.sdymo.and.abstim.lt.sdymop1) then
      jdy=l
      abstim=abstim-sdymo
c
c   Second innermost loop to find out the hour for a leap year.
c
      do k=0,23
      shrdy=k*shr
      shrdyp1=(k+1)*shr
      if(abstim.ge.shrdy.and.abstim.lt.shrdyp1) then
      jhr=k
      abstim=abstim-shrdy
c
c   Innermost loop to find out minute and second for a leap year.
c
      do n=0,59
      smnhr=n*smn
      smnhrp1=(n+1)*smn
      if(abstim.ge.smnhr.and.abstim.lt.smnhrp1) then
      jmn=n
      sec=abstim-smnhr
c
c   Date and time are setted: return.
c
      return
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      else
c
c   Second outermost loop to find out the month for a non-leap year.
c
      do m=1,12
      if(abstim.ge.smon(m).and.abstim.lt.smon(m+1)) then
      jmo=m
      abstim=abstim-smon(m)
c
c   Third outermost loop to find out the day for a non-leap year.
c
      do l=1,31
      sdymo=(l-1)*sdy
      sdymop1=l*sdy
      if(abstim.ge.sdymo.and.abstim.lt.sdymop1) then
      jdy=l
      abstim=abstim-sdymo
c
c   Second innermost loop to find out the hour for a non-leap year.
c
      do k=0,23
      shrdy=k*shr
      shrdyp1=(k+1)*shr
      if(abstim.ge.shrdy.and.abstim.lt.shrdyp1) then
      jhr=k
      abstim=abstim-shrdy
c
c   Innermost  loop to  find out  minute and second for  a  non-leap
c   year.
c
      do n=0,59
      smnhr=n*smn
      smnhrp1=(n+1)*smn
      if(abstim.ge.smnhr.and.abstim.lt.smnhrp1) then
      jmn=n
      sec=abstim-smnhr
c
c   Date and time are setted: return.
c
      return
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      endif
      enddo
c
c
      end
