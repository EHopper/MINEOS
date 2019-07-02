FFLAGS=-w -O1 -ffixed-line-length-none -fno-range-check
#LFLAGS= $(MYLFLAGS)
# FFLAGS=-ffixed-line-length-none
PROG= plot_wk
SUBS= amp.f branch_sort.f class.f color.f cvtaper.f excite.f fix_class_c.f \
	fix_class_k.f fix_class_p.f fix_class_r.f fix_class_v.f interple.f \
	interpol.f response.f search.o seek.f summary.f table.f wind.f
OBJS= $(PROG).o $(SUBS:.f=.o)

.f.o:
	gfortran $(FFLAGS) -c $*.f

#----------------------------------------------------------------------------------

$(PROG): $(OBJS) 
	gfortran $(FFLAGS) -o $(MINEOSBIN)/plot_wk $(OBJS) \
	$(MINEOSLIB)/libcip.a \
	$(MINEOSLIB)/libutil.a \
	$(MINEOSLIB)/libtau.a 

# clean up huge .o files
#	rm plot_wk.o branch_sort.o

# check object files for dependency on .h files
$(OBJS): parameter.h
	gfortran $(FFLAGS) -c $*.f
