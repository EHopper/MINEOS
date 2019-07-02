FFLAGS= $(MYFFLAGS)
LFLAGS= $(MYLFLAGS)
#
SUBS= amp.o branch_sort.o class.o color.o cvtaper.o excite_phase.o fix_class_c.o \
      fix_class_k.o fix_class_p.o fix_class_r.o fix_class_v.o interple.o \
      interpol.o response.o search.o seek.o summary.o table.o wind_phase.o
OBJS= mask_phase.o $(SUBS)


mask_phase: $(OBJS) 
	f77 $(FFLAGS) $(LFLAGS) -o $(MYBIN)/mask_phase $(OBJS) \
        $(USRLIB)/libcip.a \
        $(USRLIB)/libutil.a
# clean up huge .o files
	rm plot_wk.o branch_sort.o

