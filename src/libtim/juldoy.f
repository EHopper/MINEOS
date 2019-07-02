      integer function juldoy(iyear,idoy)
c     
c     compute Julian day from year and day of year
c     
      integer idoy, iyear, iyr, iyr1, jyear
c     
      jyear = iyear - 1900
      iyr1 = 0
      iyr = int((jyear-1)/4)
      juldoy = 2415020 + 365*jyear + idoy + iyr
c
      return
      end
