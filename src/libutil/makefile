FFLAGS=$(MYFFLAGS)
# FFLAGS= -fast -O3
#FFLAGS= $(MYFFLAGS)
LIBNAM= $(MINEOSLIB)/libutil.a
#
.f.a:
	gfortran $(FFLAGS) -c $<
	ar rv $@ $*.o
	rm -f $*.o
#
# objects for utility library
#
# 	
$(LIBNAM): \
	$(LIBNAM)(angles.o) \
	$(LIBNAM)(azimth.o) \
	$(LIBNAM)(cipget.o) \
	$(LIBNAM)(cipnum.o) \
	$(LIBNAM)(cipped.o) \
	$(LIBNAM)(cipscn.o) \
	$(LIBNAM)(cipstr.o) \
	$(LIBNAM)(cipsub.o) \
	$(LIBNAM)(ciptok.o) \
	$(LIBNAM)(ciptyp.o) \
	$(LIBNAM)(daymo.o) \
	$(LIBNAM)(doy.o) \
	$(LIBNAM)(ellip.o) \
	$(LIBNAM)(fcn.o) \
	$(LIBNAM)(gcpath.o) \
	$(LIBNAM)(gcpath_e.o) \
	$(LIBNAM)(get_bath.o) \
	$(LIBNAM)(get_unit.o) \
	$(LIBNAM)(interpol.o) \
	$(LIBNAM)(interple.o) \
	$(LIBNAM)(ival.o) \
	$(LIBNAM)(kblnk.o) \
	$(LIBNAM)(klen.o) \
	$(LIBNAM)(lpyr.o) \
	$(LIBNAM)(llen.o) \
        $(LIBNAM)(kblnk.o) \
        $(LIBNAM)(lnblnk.o) \
	$(LIBNAM)(maxsp.o) \
	$(LIBNAM)(midpnt.o) \
	$(LIBNAM)(midpnt_e.o) \
	$(LIBNAM)(moment_a.o) \
	$(LIBNAM)(moment_rtf.o) \
	$(LIBNAM)(mtimes.o) \
	$(LIBNAM)(mtimes_e.o) \
	$(LIBNAM)(pick_filter.o) \
	$(LIBNAM)(rspln.o) \
	$(LIBNAM)(sadd.o) \
	$(LIBNAM)(sdcoht.o) \
	$(LIBNAM)(sdiff.o) \
	$(LIBNAM)(sec_time.o) \
	$(LIBNAM)(sec_time_inv.o) \
	$(LIBNAM)(dsec_time.o) \
	$(LIBNAM)(dsec_time_inv.o) \
	$(LIBNAM)(sort.o) \
	$(LIBNAM)(splint.o) \
	$(LIBNAM)(spread.o) \
	$(LIBNAM)(swap.o) \
	$(LIBNAM)(ttimes.o) \
	$(LIBNAM)(xyz2geo.o)
	ranlib $(LIBNAM)
#
# dependencies
#
$(LIBNAM)(gcpath.o):		numerical.h
$(LIBNAM)(gcpath_e.o):		numerical.h
$(LIBNAM)(moment_a.o):		numerical.h
$(LIBNAM)(moment_rtf.o):	numerical.h
$(LIBNAM)(mtimes.o):		parameter.f
$(LIBNAM)(mtimes_e.o):		parameter.f
$(LIBNAM)(mtimes_e.o):		numerical.h
$(LIBNAM)(pick_filter.o):	numerical.h
$(LIBNAM)(sec_time.o):		numerical.h
$(LIBNAM)(sec_time_inv.o):	numerical.h
$(LIBNAM)(spread.o):		numerical.h
$(LIBNAM)(ttimes.o):		parameter.f
