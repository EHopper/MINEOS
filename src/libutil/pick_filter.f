      subroutine pick_filter(typbuf,ibp,ltype,w0o,w1o,w2o,w3o)
c
c     subroutine to extract info about desired filter parameters
c
      character*(*) typbuf, ltype
      character*256 sub, cmd
c
      real*4 w0o, w1o, w2o, w3o
c
      integer*4 ibp, igo
c
      include 'numerical.h'
c
      data cmd /'lobw hibw bpbw bpcos gauss cos2 butter'/
c
1005  igo = 0
      call ciptok(typbuf,ibp,sub)
      call cipscn(sub,cmd,3,igo)
      ltype = ' '
      w0o = 0.
      w1o = 0.
      w2o = 0.
      w3o = 0.
c
c     check against arg list
c
      if ((igo.eq.-3).or.(igo.eq.-1)) then
        print*,' '
        print*,' Recognized filters:'
        print*,' '
        print*,' lobw w0           - 6-pole low pass butterworth'
        print*,'                       with hi-freq cutoff w0'
        print*,' hibw w0           - 6-pole high pass butterworth'
        print*,'                       with lo-freq cutoff w0'
        print*,' bpbw w1 w2        - 6-pole band pass butterworth'
        print*,'                       between w1 and w2'
        print*,' bpcos w0 w1 w2 w3 - band pass cosine**2 filter'
        print*,'                       cosine taper between w0 and w1 '
        print*,'                   and cosine taper between w2 and w3 '
        print*,' gauss w0 w1       - gaussian filter'
        print*,'                       with center w0 and bandwidth w1'
        print*,' cos2 w0 w1        - cosine**2 filter'
        print*,'                       with center w0 and bandwidth w1'
        print*,' butter np w1 w2   - np-pole band pass butterworths'
        print*,'                       between w1 and w2'
        print*,' '
        print*,' all frequencies should be in Mhz'
        print*,' '
        go to 2000
c
c     user has wisely choosen to quit>
c
      else if ((igo.eq.-2) .or. (igo .eq. -6)) then
        go to 1099
c        
c     the program cant find a match
c
      else if (igo.eq.0) then
        print*,' i dont understand: ',token
        goto 1099
c
c     found an 'eoi'
c
      else if (igo.eq.-99) then
        goto 1099
c
c     string is too long
c
      else if(igo .eq. -999) then
        print*,' command string is too long.'
        goto 1099
c
c     go to appropriate module
c
      else
        go to (1010,1020,1030,1040,1050,1060,1070) igo
      endif
c
c     lo-pass butterworth
c
1010  continue
c
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: lobw w1 '
        call ciptyp(' enter corner frequency:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w0o = rnum
      end if
      ltype = 'lobw'
      go to 1099
c
c     hi-pass butterworth
c
1020  continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: hibw w1 '
        call ciptyp(' enter corner frequency:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w0o = rnum
      end if
      ltype = 'hibw'
      go to 1099
c
c     band pass butterworth
c
1030  continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: bpbw w1 w2'
        ibp = 0
        call ciptyp(' enter w1 and w2:',typbuf)
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w1o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter w2:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w2o = rnum
      end if
      ltype = 'bpbw'
      go to 1099
c
c     band pass cosine
c
1040  continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: bpcos w0 w1 w2 w3'
        ibp = 0  
        call ciptyp(' enter w0, w1, w2, and w3 :',typbuf) 
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w0o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter w1:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w1o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter w2:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w2o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter w3:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w3o = rnum
      end if
      ltype = 'bpcos'
      go to 1099
c
c     gaussian filter
c
1050  continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: gauss wo sigma'
        ibp = 0
        call ciptyp(' enter wo and sigma:',typbuf)
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w0o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter sigma:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w1o = rnum
      end if
      ltype = 'gaus'
      go to 1099
c
c     cosine**2 filter
c
1060  continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: cos2 wo sigma'
        ibp = 0
        call ciptyp(' enter wo and sigma:',typbuf)
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w0o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter sigma:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w1o = rnum
      end if
      ltype = 'cos2'
      go to 1099
c
c     n-pole band pass butterworth
c
1070  continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        print*,' command format: butt np w1 w2'
        ibp = 0
        call ciptyp(' enter np and w1 and w2:',typbuf)
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag) 
      if (iflag.eq.0) then 
        npole = rnum 
      end if 
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        ibp = 0
        call ciptyp(' enter w1 and w2:',typbuf)
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w1o = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter w2:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        w2o = rnum
      end if
      ltype = 'butt'
      go to 1099
c
c     all done
c
1099  continue
c
c     check for different units
c
      call ciptok(typbuf,ibp,sub)
      if (sub .ne. 'eoi') then
        if (sub .eq. 's') then
          if (w0o .ne. 0.0) then
            w0o = 1000.0 / w0o
          endif
          if (w1o .ne. 0.0) then
            w1o = 1000.0 / w1o
          endif
          if (w2o .ne. 0.0) then
            w2o = 1000.0 / w2o
          endif
          if (w3o .ne. 0.0) then
            w3o = 1000.0 / w3o
          endif
          if (ltype.eq.'gaus' .or. ltype.eq.'cos2') then
            chs = (1000./w0o + 1000./w1o)
            cls = (1000./w0o - 1000./w1o)
            clf = 1000./cls
            chf = 1000./chs
            w1o =  clf - chf
            print*,'Guass params in mHz: ',w0o,w1o
          endif
        elseif (sub .eq. 'h') then
          if (w0o .ne. 0.0) then
            w0o = 1000.0 * w0o
          endif 
          if (w1o .ne. 0.0) then
            w1o = 1000.0 * w1o
          endif
          if (w2o .ne. 0.0) then
            w2o = 1000.0 * w2o
          endif
          if (w3o .ne. 0.0) then
            w3o = 1000.0 * w3o
          endif
        elseif (sub .eq. 'r') then
          if (w0o .ne. 0.0) then
            w0o = rad * w0o
          endif 
          if (w1o .ne. 0.0) then
            w1o = rad * w1o
          endif
          if (w2o .ne. 0.0) then
            w2o = rad * w2o
          endif
          if (w3o .ne. 0.0) then
            w3o = rad * w3o
          endif
        endif
      endif
      if (ltype .eq. 'butt') then
        w0o = npole
      end if
2000  continue
      return
      end
