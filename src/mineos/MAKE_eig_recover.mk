FFLAGS=-w -O1 -ffixed-line-length-none -fno-range-check
#LFLAGS= $(MYLFLAGS)
LFLAGS=
PROG= eig_recover
SUBS= kblnk.f
OBJS= $(PROG).o $(SUBS:.f=.o)


.f.o:
	gfortran $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS)
	gfortran $(FFLAGS) $(LFLAGS) -o $(MINEOSBIN)/$@ $(OBJS)

# check object files for dependency on .h files
$(OBJS): parameter.h
	gfortran $(FFLAGS) -c $*.f
