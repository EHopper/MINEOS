      double precision function tcalc(iy,id,ih,im,sec)
c     
c     calculate time
c     
      real*8 jtime
      real*4 sec
      integer jdcalc, juldoy, iy, id, ih, im
      jdcalc = juldoy(iy,id)
      tcalc  = jtime(jdcalc,ih,im,sec)
c
      return
      end
