#  Compiler options.
FC = gfortran
FFLAGS = $(MYFFLAGS) 
#For HP
# FFLAGS = $(MYFFLAGS) +O3 +U77 +E5 
#
#  Library pathname.
#
LIBNAM= $(MINEOSLIB)/libtim.a
#
#  Compile, archive and clean.
#
.f.a:
	$(FC) $(FFLAGS) -c  $<
	ar rv $@ $*.o
	rm -f $*.o
#
#  List all the target objects.
#
$(LIBNAM): \
	$(LIBNAM)(abstime.o) \
	$(LIBNAM)(jday.o) \
	$(LIBNAM)(jtime.o) \
	$(LIBNAM)(juldoy.o) \
	$(LIBNAM)(tcalc.o) \
	$(LIBNAM)(time.o) \
	$(LIBNAM)(timeph.o) \
	$(LIBNAM)(timej.o)
#
#  Set index.
#
$(LIBNAM): ; ranlib $(LIBNAM)
