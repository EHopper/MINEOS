FFLAGS= $(MYFFLAGS)
LIBNAM= $(MINEOSLIB)/libcip.a
.f.a:
	gfortran $(FFLAGS) -c  $<
	ar rv $@ $*.o
	rm -f $*.o
#
#  List all the target objects
#
# $(LIBNAM): \
# 	$(LIBNAM)(ciplib.o)
# 		ranlib $(LIBNAM)
#
$(LIBNAM): \
	$(LIBNAM)(ciplib.o) 
#
#  Set index.
#
$(LIBNAM): ; ranlib $(LIBNAM)
