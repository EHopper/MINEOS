c
c
      subroutine timej(abstime,jyr,jjd,jhr,jmn,sec)
c
c   Purpose:
c     To transform  a period of time, expressed as seconds from 1900
c     1 0 0 0.0, in a date (julian day) and time.     
c
c   WARNINGS: - the routine works only up to the year 1999.
c             - the routine  does not have a control  on the century
c               and for istance it assumes 1979 to be 79.
c
c   02/02/03 -- modified to return 4-digit year -- above two warnings
c               obsolete.  see comments for syr in catbin program
c               (usr/local/src/Cmt).  Currently up to date through 2004.
c   Arguments:
c     abstime    time in seconds from 1900 1 0 0 0.0.
c     jyr        year.
c     jjd        day of the year (julian day).
c     jhr        hour.
c     jmn        minute.
c     sec        second.
c
c
      implicit real*8 (a-h,o-z)
      dimension syr(0:110)
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
      data (syr(i),i=49,110)/
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
      data sjd/0.86400d+05/
      data shr/0.3600d+04/
      data smn/0.60d+02/
c
      abstim=abstime
c
c   Outermost loop to find out the year.
c
      do i=0,105
      if (abstim.ge.syr(i).and.abstim.lt.syr(i+1)) then
      jyr=1900+i
      abstim=abstim-syr(i)
c
c   Second outermost loop to find out the day of the year.
c
      do l=1,366
      sjdyr=(l-1)*sjd
      sjdyrp1=l*sjd
      if(abstim.ge.sjdyr.and.abstim.lt.sjdyrp1) then
      jjd=l
      abstim=abstim-sjdyr
c
c   Second innermost loop to find out the hour.
c
      do k=0,23
      shrjd=k*shr
      shrjdp1=(k+1)*shr
      if(abstim.ge.shrjd.and.abstim.lt.shrjdp1) then
      jhr=k
      abstim=abstim-shrjd
c
c   Innermost loop to find out minute and second.
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
c
c
      end
