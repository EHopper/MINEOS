#FFLAGS= $(MYFFLAGS)
#LFLAGS= $(MYLFLAGS)
FFLAGS=-w -O1 -ffixed-line-length-none -fno-range-check
#LFLAGS= -Bstatic
LFLAGS=
#
PROG= mineos_nohang
SUBS= baylis.f bfs.f dermf.f derms.f detqn_nohang.f drspln.f dsplin.f eifout.f entry.f\
      fprop.f fprpmn.f fpsm.f fsbdry.f fsbm.f gauslv.f grav.f intgds.f match.f\
      model.f modout.f ortho.f remedy_nohang.f rkdot.f rotspl_nohang.f rprop.f rps.f sdepth.f\
      sfbdry.f sfbm.f sprop.f sprpmn.f spsm.f startl.f steps.f svd.f tprop.f\
      tps.f trknt.f whead.f wtable.f zknt.f
OBJS= mineos.o $(SUBS:.f=.o)
# FC=gfortran

.f.o:
	gfortran $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	gfortran $(FFLAGS) $(LFLAGS) -o $(MINEOSBIN)/$@ $(OBJS)

# check object files for dependency on .h files
$(OBJS): parameter.h
	gfortran $(FFLAGS) -c $*.f

clean:
	-rm *.o 