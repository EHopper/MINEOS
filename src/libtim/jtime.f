      double precision function jtime(jday,jhr,jmin,sec)
c     
c     calculate time in days from jday,jhr,jmin,jsec
c     
      real*8 tempo
      integer jday,jhr,jmin
      real sec
      jtime = dble(jday) 
     &     + ((((dble(sec)/60.d0)+dble(jmin))/60.d0)+dble(jhr))/24.d0
c
      return
      end
