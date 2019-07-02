OBJTAU = libtau.o libsun.o
FFLAGS=$(MYFFLAGS)
#FFLAGS = -O3
LFLAGS = -L$(MINEOSLIB)
LLIBS = -ltau
ARFLAGS= rcv

libtau.a:	${OBJTAU}
		ar ${ARFLAGS} $(MINEOSLIB)/libtau.a ${OBJTAU}
		ranlib $(MINEOSLIB)/libtau.a
setbrn:		setbrn.o
		gfortran $(FFLAGS) $(MINEOSLIB)/libtau.a setbrn.o -o $(MINEOSBIN)/setbrn
ttimes:		ttimes.o
		gfortran $(FFLAGS) $(LFLAGS) ttimes.o $(LLIBS) -o $(MINEOSBIN)/ttimes


$(OBJTAU): libtau.f libsun.f
	gfortran $(FFLAGS) -c $*.f
