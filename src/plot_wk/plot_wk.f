      program plot_wk
c
c     program to plot wk diagram using gpr routines - and much much more
c     based on command processor - the main program is rather long...
c     uses gpr plotting calls
c 
c     written November 1988 by Lind Gee
c     06/96/89 changed to make fill matrix real, allowing for amplitude windowing
c     07/18/90 plotting of frechet kernels added
c     08/03/90 plotting of eigenfunctions added  
c     02/19/92 Sun version - no plotting
c     05/23/02 Add phname optoin == choose group-vel window
c               by phase name.  Requires new version of exc
c               files that include source depth in header
c
      include 'parameter.h'
      include 'parameter.f'
c      parameter (nph=10)
c
      real*4 modes(0:maxmodes,nprop)
      real*4 wmin, wmax, pmin, pmax, gmin, gmax
      real*4 xin(2), xout(2), xtic
      real*4 yin(2), yout(2), ytic
      real*4 xlow, xhigh, ylow, yhigh
      real*4 pvel(10)
      real*4 fill(0:maxmodes)
      real*4 fmin, fmid, fmax
      real*4 stuff(4)
      real*4 c0, c1, c2, c3
      real*4 dist,dep,tt(nphase)
      real*4 dtdd(nphase),dtdh(nphase),dddp(nphase),phcd(nphase)
c
      integer*4 llmin, llmax
      integer*4 nnmin, nnmax
      integer*4 llow, lup, nlow, nup
      integer*4 wtemp, index, jindex
      integer*4 xlen, ylen
      integer*4 s1(100), s2(100)
      integer*4 nnb(0:maxl,nbranch), llb(0:maxl,nbranch), numb(nbranch)
      integer*4 iprop(4)
c
      character*1 ans, comp(3), ch1, ch, wtype
      character*2 ch2
      character*3 ch3a, ch3b
      character*4 ch4
      character*5 ch5
      character*6 ch6
      character*8 pname
      character*10 cht
      character*80 a_desc
      character*256 text, text1, text2, text3, text4, text5
      character*256 xlabel, ylabel, title
      character*256 typbuf, token, sub                                                            
      character*256 cmd, cmd2
      character*256 m_file,o_file,c_file,w_file,b_file,e_file,r_file
      character*256 f_file
c
c     gpr declarations
c
      integer*2 config, unit
      integer*2 disp(31), disp_len, disp_len_ret
      integer*2 bmwdth, bmhght
      integer*2 init_bitmap_size(2), hi_plane_id
      integer*2 disp_bitmap_size(2)
      integer*2 xoff, yoff
      integer*2 x, y, x1, y1
      integer*2 box(4)
      integer*2 pause(3)
      integer*2 center(2), radius
      integer*2 setsize, ev_pos(2), ev_type, keys(40)
      integer*2 position(2)
      integer*2 cursor_op
      integer*2 font_id_sm, font_id_old, direction, font_id_cur
      integer*2 clip(4)
      integer*2 name_size, version(2), groups, group_header(8)
      integer*2 window(2,2), dest(2)
      integer*2 size(2)
      integer*2 ixlim, iylim, ixoff, iyoff
c
      integer*4 origin(2)
      integer*4 cursor_pat
      integer*4 status
      integer*4 init_bitmap
      integer*4 c_value(1:16), icol(0:15)
      integer*4 inx(0:maxl,0:maxn)
      integer*4 ifill(0:maxmodes),cfill(0:maxmodes)
      integer*4 icolor_fill
      integer*4 inx2(0:maxmodes,2)
      integer*4 icolor, ncolor, mcolor
cad      integer*4 gpr_$light_blue, gpr_$light_purple
      integer*4 attrib_block_desc, save_bitmap
c
      character ev_char, bitmap_name*80
      logical in_color, wait, active, line, lhelp, lcreate
      logical lhigh, lsummary, lpv, lselect, lread
      logical lomega, lphase, lgroup, lfirst, lcolor, lbranch
      logical lkernel, lwin, lphn
c
      common /bitmap/ bmwdth, bmhght, xoff, yoff
      common /branch/ lbranch, nbr, numb, nnb, llb
      common /lost/ lost(0:maxl,0:maxn)
      common /c_color/ icol
      common /color2/ fmin, fmid, fmax
      common /color3/ lcolor, clmin, clmax
      common /fill/ lsummary, fill, ifill, cfill
      common /limits/ jcom, nmodes, nlow, nup, llow, lup
      common /limits2/ wmin, wmax, pmin, pmax, gmin, gmax
      common /mode/ inx, modes, inx2
      common /nparam/ nparam, nextra
      common /plot/ lomega,lphase,lgroup,xlow,xhigh,ylow,yhigh,
     &              xlen,ylen,ix,iy 
      common /plot2/ xlabel, ylabel
      common /c_search/ isearch, index(nprop), value(nprop,2)
      common /c_excite/ dist, dep
c
c      data xoff,yoff / 100, 100/
      data ixoff,iyoff / 100, 100/
      data position / 512,500/
      data pause / 0,20,0/
      data tpi /6.2831853071796/
      data comp /'S','T','S'/
      data lhelp, lhigh / 2*.false. /
      data version /1,1/
      data groups / 1/
      data group_header / 4, 1, 0, 0, 0, 0, 0, 0/
      data bitmap_name, name_size /'wk_bm', 5/
      data ibm / 0/
      data title /' '/
      data lpv /.false./
      data pvel /4.53,6.58,8.23,11.94,13.26,14.76,18.07,25.08,57.37,0./
      data npvel /9/
c      data cfill /maxmodes1*2/
c      data ifill /maxmodes1*0/
      data mcolor /16/
c      data nparam /4/
      data lread /.false./
c      data lomega, lphase, lgroup /.true.,.false.,.false./
c      data ylen, ylabel /15, 'Frequency (mHz)'/
c      data xlen, xlabel /13, 'Angular Order'/
c      data ix, iy /1, 1/
      data lfirst /.true./
      data icolor_fill /8/
c      data nextra /0/
      data lkernel /.false./
      data lwin,lphn /2*.false./
c
      data cmd /'summary table search output clear ctable title phase wi
     &ndow fill old branch excite type resp class color fre subt eig cvt
     &aper boff outbran phname plot stop'/

      include 'numerical.h'
	  
c     initialize variables in common blocks
      xoff=100
	  yoff=100
	  lsummary=.false.
	  nparam=4
	  lomega=.true.
	  lgroup=.false.
	  lphase=.false.
	  ylen=15
	  ylabel='Frequency (mHz)'
	  xlen=13
	  xlabel='Angular Order'
	  ix=1
	  iy=1
	  nextra=0
c
c     establish initial color table
c  
c     c_value(ii) = ir*65536 + ig*256 + ib
c
cad     c_value(1) = gpr_$black
cad     c_value(2) = gpr_$white
c
c     mode = 'off'
c
cad      c_value(3) = gpr_$blue
c   
c     table for amplitude
c
c     old table
c
c      c_value(4) =              178*256 + 255
c      c_value(5) =  255*65536
c      c_value(6) =  255*65536 + 120*256  
c      c_value(7) =  255*65536 + 160*256 
c      c_value(8) =  255*65536 + 200*256 
c
c     new table  - only values 5, 6, 8, and 9 are used
c
      c_value(4) =              178*256 + 255
      c_value(5) =              200*256 + 255
      c_value(6) =  255*65536
      c_value(7) =  255*65536 + 160*256 
      c_value(8) =  255*65536 + 160*256 
c
c     mode = 'on'
c
      c_value(9) = 255*65536 + 255*256 
c
c     mode = highlighted'
c
cad      c_value(10) = gpr_$green
c
c     white
c
      do ii = 11, 16
cad        c_value(ii) = gpr_$white
      end do
c
      radius = 1
      unit = 1
      scale = 0.5
c
c     main loop of command processor
c
9000  continue
      call ciptyp('command:',typbuf)
      ibp = 0
9010  igo = 0
      call ciptok(typbuf,ibp,token)
      call cipscn(token,cmd,3,igo)
c
c  check against arg list
c
c     print list of recognized commands
c
      if ((igo .eq. -3) .or. (igo .eq. -1)) then
        print*,' '
        print*,' Recognized commands:'
        print*,'   summary - read mode summary file'
        print*,'   table   - read mode table header file'
        print*,'   search  - set search table'
        print*,'   output  - specify name of outputfile for fill matrix'
        print*,'   clear   - clear fill table and reset parameters'
c       print*,'   ctable  - specify color table'
c       print*,'   title   - set title for plot'
c       print*,'   phase   - set new limits for constant phase velocity'
        print*,'   window  - calculate group velocity window (gaussian)'
        print*,'   phname  - calculate group velocity window (gaussian)'
        print*,'             for a specified phase (P, S, SS, etc). '
        print*,'             Needs width, corr (phname S 0.01 0) '
        print*,'   fill    - read a pre-existing fill table'
        print*,'   old     - specify old summary file'
        print*,'   branch  - read a pre-existing branch file'
        print*,'   excite  - read an excitation file'
c       print*,'   type    - specify plot type: w, c, or g'
        print*,'   resp    - add instr group delay for windowing'
        print*,'   class   - classify spheroidal modes: K,C,V,R'
c       print*,'   color   - specify range of values for color'
        print*,'   fre     - specify file for frechet kernels'
        print*,'             eig and fre are exclusive options'
        print*,'   subt    - subtract selected modes from mode set'
        print*,'   eig     - specify file for eigenfunctions'
        print*,'             eig and fre are exclusive options'
        print*,'   cvtaper - taper in phase velocity'
        print*,'             cos**2 taper between c0 & c1 and c2 & c3'
        print*,'             specify 4 values in km/s'
        print*,'   boff    - branch off - turn off particular branches'
        print*,'   outbran - output multiple-branch masks'
        print*,'             Enter first and last branch in range'
c       print*,'   plot - plot'
        go to 9000
c
c     quit this program
c
      else if ((igo .eq. -2) .or. (igo .eq. -6)) then
        go to 9999
c
c     unrecognized command
c
      else if (igo .eq. 0) then
        print*,' i dont understand: ',token
        go to 9010
c
c     found an 'eoi'
c
      else if (igo .eq. -99) then
        goto 9000
c
c     string is too long
c     
      else if(igo .eq. -999) then
        print*,' command string is too long.'
        goto 9010
c
c     go to appropriate module
c
      else
        go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     &         170,180,190,200,210,220,230,240,999,9999) igo
      endif
c
c     entry point for command summary
c
10    continue
      if (nparam .eq. 4) then
        nparam = 25
      endif
      call ciptok(typbuf,ibp,m_file)
      if (m_file .eq. 'eoi') then
        call ciptyp(' enter name of mode summary:',m_file)
      endif
      call summary(m_file)
      nnmax = nup
      nnmin = nlow
      llmax = lup
      llmin = llow
      wwmin = wmin
      wwmax = wmax
      do ii = 1, nmodes
        fill(ii) = 0.
        ifill(ii) = 0
      end do
      lselect = .false.
      lsummary = .true.
      if (lomega) then
        xlow   = real(llow)
        xhigh  = real(lup)
        ylow   = wmin
        yhigh  = wmax
      elseif (lphase) then
        xlow   = pmin
        xhigh  = pmax
        ylow   = wmin
        yhigh  = wmax
      elseif (lgroup) then
        lgroup = .true.
        lomega = .false. 
        lphase = .false. 
        xlow   = gmin
        xhigh  = gmax
        ylow   = wmin
        yhigh  = wmax
      endif      
c
      go to 9010
c
c     entry point for command table
c       note the difference between the commands table and summary
c
20    continue
      if (nparam .gt. 4) then
        nparam = 4
      endif
      call ciptok(typbuf,ibp,m_file)
      if (m_file .eq. 'eoi') then
        call ciptyp(' enter name of mode table header:',m_file)
      endif
      call table(m_file)
      nnmax = nup
      nnmin = nlow
      llmax = lup
      llmin = llow
      wwmin = wmin
      wwmax = wmax
      do ii = 1, nmodes
        fill(ii) = 0.
        ifill(ii) = 0
      end do
      lselect = .false.
      lsummary = .false.
      if (lomega) then
        xlow   = real(llow)
        xhigh  = real(lup)
        ylow   = wmin
        yhigh  = wmax
      elseif (lphase) then
        xlow   = pmin
        xhigh  = pmax
        ylow   = wmin
        yhigh  = wmax
      elseif (lgroup) then
        lgroup = .true.
        lomega = .false. 
        lphase = .false. 
        xlow   = gmin
        xhigh  = gmax
        ylow   = wmin
        yhigh  = wmax
      endif      
c
      go to 9010
c
c     entry point for command search
c
30    continue
      call search(nm,lselect)
      print*, nm, ' modes selected by search'
      go to 9010
c
c     entry point for command output
c
40    continue
      call ciptok(typbuf,ibp,o_file)
      if (o_file .eq. 'eoi') then
        call ciptyp(' enter name of output file:',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,o_file)
      endif
      call ciptok(typbuf,ibp,sub)
c
c     the 'b' option allows a specified branch to be output
c     will overwrite fill matrix
c
      if (sub(1:1) .eq. 'b') then
        call ciptyp(' specifiy branch number (0=Fundamental): ',sub)
        call cipnum (sub, rnum, iflag)   
        ib = rnum + 1
        do jj = 1, nmodes
          fill(jj) = 0.0
          ifill(jj) = 0
          cfill(jj) = 0
        end do
        do jj = 1, numb(ib)
          nn = nnb(jj,ib)
          ll = llb(jj,ib)
          ind = inx(ll,nn)
          fill(ind) = 1.0
          ifill(ind) = 1
        end do
        call color(nmodes)
      endif
c
c     open binary output file and write out fill matrix
c
      open(unit=4,file=o_file,form='unformatted',access='sequential')
      write(4) jcom, nmodes, nlow, nup, llow, lup, wgrav
      print*,' enter ascii descriptor (<80 char)'
      read(*,'(a80)') a_desc
      print*,' reading fill table for: ',nmodes,' modes'
      print*, a_desc
      print*,' search algorithm: '
      do ii = 1, isearch
        print*,ii,index(ii),value(index(ii),1),value(index(ii),2)
      end do
      if (gv0 .ne. 0.) then
        print*,' window algorithm: ',gv0, s0
      endif
      write(4) a_desc
      write(4) isearch, (index(ii),(value(index(ii),jj),jj=1,2),
     +         ii=1,isearch)
      write(4) gv0, s0
      do ii = nlow, nup
        write(4) ii,(fill(inx(jj,ii)), jj = llow, lup)
      end do
      close(4)
      go to 9010
c
c     entry point for command clear
c
50    continue
c
c     clear selected modes by setting fill = 0.
c
      do ii = 0, nmodes
        fill(ii) = 0.
        ifill(ii) = 0
        cfill(ii) = 2
      end do
      lselect = .false.
      do ii = 1, nprop
        index(ii) = 0
        value(ii,1) = 0.0
        value(ii,2) = 0.0
      end do
      isearch = 0
      s0 = 0.0
      gv0 = 0.0
      fmin = 0.0
      fmid = 0.0
      fmax = 0.0
      clmin = 0.0
      clmax = 0.0
      lcolor = .false.
      title = ' '
      ititle = 0
      scale = 0.5
      go to 9010
c
c     entry point for command ctable
c
60    continue
      call ciptok(typbuf,ibp,c_file)
      if (c_file .eq. 'eoi') then
        call ciptyp(' enter name of color table:',c_file)
      endif
      open(4,file=c_file)
      read(4,*) ncolor
      do ii = 1, ncolor
        read(4,*) ir, ig, ib
        c_value(ii) = ir*65536 + ig*256 + ib
      end do
      close(4)
      go to 9010      
c
c     entry point for command title
c
70    continue
      call ciptyp(' enter title of plot:',title)
      ii = 256
      do while (title(ii:ii) .eq. ' ')
        ii = ii - 1
      end do
      ititle = ii
      go to 9000
c
c     entry point for command phase
c
80    continue
      print*,'  #  Phase velocity'
      do ii = 1, npvel
        print*, ii, pvel(ii)
      end do
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter # of lines and corresponding velocities:',
     &    typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub, rnum, iflag)   
      npvel = rnum
      if (npvel .gt. 10) then
        npvel = 10
      end if
      do ip = 1, npvel
        call ciptok(typbuf,ibp,sub)
        call cipnum (sub, rnum, iflag)   
        pvel(ip) = rnum
      end do
      print*,' '
      print*,'  #  Phase velocity'
      do ii = 1, npvel
        print*, ii, pvel(ii)
      end do
      lpv = .true.
      go to 9000
c
c     entry point for window
c
90    continue
      if (lphn) then
        print*,'Already windowed using phname -- wind ignored'
        print*,''
        go to 9000
      endif
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter group velocity for window peak (km/sec):',
     &              sub)
      endif
      call cipnum (sub, rnum, iflag)   
      gv0 = rnum
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter window width (1/sec):',sub)
      endif
      call cipnum (sub, rnum, iflag)   
      s0 = rnum
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter parameters in group delay (0=group; 1=source
     &; 2=instr):',sub)
      endif
      call cipnum (sub, rnum, iflag)   
      np = rnum
      print*,'gv0,s0,np',gv0,s0,np
      call wind(gv0, s0, np)
      lwin = .true.
      go to 9000
c
c     entry point for command fill
c
100   continue
      call ciptok(typbuf,ibp,o_file)
      if (o_file .eq. 'eoi') then
        call ciptyp(' enter name of fill file:',o_file)
      endif
c
c     open binary output file and read in fill matrix
c
      open(unit=4,file=o_file,form='unformatted',access='sequential',
     +     status='old')
      read(4) jcom2, nmodes2
      read(4) a_desc
      read(4) isearch, (index(ii),(value(ii,jj),jj=1,2),ii=1,isearch)
      read(4) gv0, s0
      print*,' reading fill table for: ',nmodes2,' modes'
      print*, a_desc
      print*,' search algorithm: '
      do ii = 1, isearch
        print*,ii,index(ii),value(ii,1),value(ii,2)
      end do
      if (gv0 .ne. 0.) then
        print*,' window algorithm: ',gv0, s0
      endif
      do ii = nlow, nup
        read(4) idum,(fill(inx(jj,ii)), jj = llow, lup)
        do jj = llow, lup
          if (fill(inx(jj,ii)) .ne. 0.0) then
            ifill(inx(jj,ii)) = 1
            print*,' filling mode: ', ii, jj, fill(inx(jj,ii))
          else
            ifill(inx(jj,ii)) = 0
          end if
        end do
      end do
      call color(nmodes)
      close(4)
      go to 9010
c
c     entry point for command old (toggles nparam value)
c
110   continue
      if (nparam .gt. 13) then
        nparam = 13
      else
        nparam = 25
      endif
      go to 9010
c
c     entry point for command branch
c
120   continue
      call ciptok(typbuf,ibp,b_file)
      if (b_file .ne. 'eoi') then
        lread = .true.
      else
c
c       give them another chance
c 
        call ciptyp(' enter branch file or <ret> to sort:',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,b_file)
        if (b_file(1:3) .ne. 'eoi') then
          print*,' be sure to have selected all the desired modes'
          lread = .true.
        endif
      endif
c
c     read a pre-existing file
c
      if (lread) then
        if (b_file(1:5) .eq. 'clear') then
          do ii = 1, nbranch
            do jj = 0, maxl
              nnb(jj,ii) = 0
              llb(jj,ii) = 0
            end do
          end do
          nbr = 0
          do ii = 1, nmodes
            fill(ii) = 0.
            ifill(ii) = 0
            cfill(ii) = 2
          end do
        else
          open(4,file=b_file,form='unformatted',access='sequential')
          read(4) 
          read(4) nbr2, (numb(nbr + ii), ii = 1, nbr2)
          do ii = 1, nbr2
            i2 = nbr + ii
            read(4) (nnb(jj,i2), llb(jj,i2), 
     &              (stuff(kk),kk = 1,4),jj = 1,numb(i2))
            do jj = 1, numb(i2)
              nn = nnb(jj,i2)
              ll = llb(jj,i2) 
              ind = inx(ll,nn)
              if (ind .eq. 0) then
                print*,' zero index: ', nn, ll, ii, jj
              endif
            end do   
          end do
          close(4)
          nbr = nbr + nbr2
          print *,' '
          print *,' information for ',nbr2,' branches read'
          print *,' total of ',nbr,' branches stored'
          lbranch = .true.
        endif  
c
c      sort modes selected by the 'class' commands
c
      elseif (.not.(lread)) then
        call kblnk(m_file,kk)
        b_file = m_file(1:kk)//'.branch'
        if ((jcom .eq. 2) .or. (icolor_fill .eq. 8)) then
ccc     order by n
          call branch_sort(b_file)
        elseif (icolor_fill .eq. 4) then
ccc     k modes
          call fix_class_k (b_file,icolor_fill)
        elseif (icolor_fill .eq. 6) then
ccc     r modes
          call fix_class_r (b_file,icolor_fill)
        elseif (icolor_fill .eq. 3) then
ccc     c modes
          call fix_class_c (b_file,icolor_fill)
        elseif (icolor_fill .eq. 5) then
ccc     v modes
          call fix_class_v (b_file,icolor_fill)
        elseif (icolor_fill .eq. 7) then
ccc     p modes
          call fix_class_p (b_file,icolor_fill)
        endif
      endif
      lread = .false.
      go to 9010
c
c     entry point for command excite
c
130   continue
      call ciptok(typbuf,ibp,e_file)
      if (e_file .eq. 'eoi') then
        call ciptyp(' enter name of excitation file:',e_file)
      endif
      call excite(e_file,lread)
      go to 9010
c
c     entry point for command type
c
140   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter plot type: w, c, g',sub)
      endif
c     call plot_type(sub(1:1))
      lfirst = .true.
      go to 9010
c
c     entry point for command resp
c
150   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter name of response file: ',sub)
      endif
      r_file = sub
      call response(r_file)
      go to 9010
c
c     entry point for command class
c
160   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter class type: ',sub)
      endif
      call class(sub(1:1),icolor_fill)
      lselect = .true.
      go to 9010
c
c     entry point for command color
c
170   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter min and max values for color: ',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub, rnum, iflag)   
      clmin = rnum
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter max value for color: ',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub, rnum, iflag)   
      clmax = rnum
      lcolor = .true.      
      call color(nmodes)
      go to 9010
c
c     entry point for command fre
c
180   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter name of frechet file: ',sub)
      endif
      f_file = sub
      lkernel = .true.
      lfre = 1
      go to 9010
c
c     entry point for command subt
c
190   continue
      do ii = 1, nmodes
        if (cfill(ii) .gt. 2) then
          fill(ii) = 0.
          ifill(ii) = 0
          cfill(ii) = 2
        else
          cfill(ii) = 8
          ifill(ii) = 1
          fill(ii) = 1.
        endif
      end do
      lselect = .true.
      go to 9010
c
c     entry point for command eig
c
200   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter name of eigenfunction file: ',sub)
      endif
      f_file = sub
      lkernel = .true.
      lfre = 2
      go to 9010
c
c     entry point for command cvtaper
c
210   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' command format: cvtaper c0 c1 c2 c3',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag .eq. 0) then
        c0 = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter c1, c2, and c3:',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag .eq. 0) then
        c1 = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter c2 and c3:',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        c2 = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter c3:',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        c3 = rnum
      end if
      call cvtaper(c0,c1,c2,c3)
      go to 9010
c
c     entry point for command boff
c
220   continue
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' command format: boff n0 n1',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag .eq. 0) then
        n0 = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter n1',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag .eq. 0) then
        n1 = rnum
      end if
      if ((n1+1) .gt. nbranch) then
        n1 = nbranch - 1
      endif
      do kk = n0 + 1, n1 + 1
        do jj = 1, numb(kk)
          nn = nnb(jj,kk)
          ll = llb(jj,kk)                
          ind = inx(ll,nn)
          fill(ind) = 0.0
          ifill(ind) = 0
        end do
      end do
      go to 9010
c
c     entry point for command outbran
c
230   continue
      call ciptok(typbuf,ibp,o_file)
      if (o_file .eq. 'eoi') then
        call ciptyp(' enter name of output file:',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,o_file)
      endif
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' specify 1st branch (0=Fundamental): :',typbuf)
        ibp = 0
        call ciptok(typbuf,ibp,sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        nob1 = rnum
      end if
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' specify last branch: :',sub)
      endif
      call cipnum (sub,rnum,iflag)
      if (iflag.eq.0) then
        nob2 = rnum
      end if
c
c     output selected branches. 
c     will overwrite fill matrix
c
      if (nob1.ge.0 .and. nob2.le.(nbr-1) .and. nob1.le.nob2) then
        do jj = 1, nmodes
          fill(jj) = 0.0
          ifill(jj) = 0
          cfill(jj) = 0
        end do
        do ib = nob1+1, nob2+1
          do jj = 1, numb(ib)
            nn = nnb(jj,ib)
            ll = llb(jj,ib)
            ind = inx(ll,nn)
            fill(ind) = 1.0
            ifill(ind) = 1
          end do
        end do
        call color(nmodes)
      endif
c
c     open binary output file and write out fill matrix
c
      open(unit=4,file=o_file,form='unformatted',access='sequential')
      write(4) jcom, nmodes, nlow, nup, llow, lup, wgrav
      print*,' enter ascii descriptor (<80 char)'
      read(*,'(a80)') a_desc
      print*,' reading fill table for: ',nmodes,' modes'
      print*, a_desc
      print*,' search algorithm: '
      do ii = 1, isearch
        print*,ii,index(ii),value(index(ii),1),value(index(ii),2)
      end do
      if (gv0 .ne. 0.) then
        print*,' window algorithm: ',gv0, s0
      endif
      write(4) a_desc
      write(4) isearch, (index(ii),(value(index(ii),jj),jj=1,2),
     +         ii=1,isearch)
      write(4) gv0, s0
      do ii = nlow, nup
        write(4) ii,(fill(inx(jj,ii)), jj = llow, lup)
      end do
      close(4)
      go to 9010
c
c     entry point for phname
c
240   continue
      if (lwin) then
        print*,'Already windowed using wind -- phname ignored'
        print*,''
        go to 9000
      endif
      call ciptok(typbuf,ibp,pname)
      if (pname .eq. 'eoi') then
        call ciptyp(' enter phase name for window peak (e.g. P,S,SS):',
     &              pname)
      endif
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter window width (1/sec):',sub)
      endif
      call cipnum (sub, rnum, iflag)   
      s0 = rnum
      call ciptok(typbuf,ibp,sub)
      if (sub .eq. 'eoi') then
        call ciptyp(' enter parameters in group delay (0=group; 1=source
     &; 2=instr):',sub)
      endif
      call cipnum (sub, rnum, iflag)   
      np = rnum

      dist=dist/dkm
      call ttimes(dep,dist,pname,nnp,tt,dtdd,dtdh,dddp,phcd)
      dist=dist*dkm

c     apply arbitrary time shift to get closer to peak of arrival

      gv0 = dist/(tt(1)+10.)
      print*,'GV Window based on phase: ',pname,' ',gv0,s0,np
c      print*, pname,dist,dep, tt(1)
      call wind(gv0, s0, np)
      lphn = .true.
      go to 9000
c
c     entry point for command plot
c
999   continue
      go to 9010
c
c
c 
c
9999  continue
      if (lkernel) then
        close(4)
      endif
      stop
      end
c
c
c
