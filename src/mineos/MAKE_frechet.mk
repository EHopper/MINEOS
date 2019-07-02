#LFLAGS= -L$(MINEOSLIB)
# FFLAGS=-ffixed-line-length-none
FFLAGS=-w -O1 -ffixed-line-length-none -fno-range-check
LDLIBS=-lgfortran -lm

frechetcomp = draw_frechet_gv frechet frechet_gv frechet_cv draw_frechet_gv

#----------------------------------
#----------------------------------

all:  clean $(frechetcomp)

#----------------------------------
#----------------------------------


frechet: $(MINEOSBIN)/frechet

$(MINEOSBIN)/frechet: frechet.o
	-rm $(MINEOSBIN)/frechet
	gfortran $(FFLAGS) -o $(MINEOSBIN)/frechet frechet.o

frechet.o: frechet.f
	gfortran $(FFLAGS) -c -o frechet.o frechet.f

#----------------------------------

frechet_gv: $(MINEOSBIN)/frechet_gv

$(MINEOSBIN)/frechet_gv: frechet_gv.o
	-rm $(MINEOSBIN)/frechet_gv
	gfortran $(FFLAGS) -o $(MINEOSBIN)/frechet_gv frechet_gv.o

frechet_gv.o: frechet_gv.f
	gfortran $(FFLAGS) -c -o frechet_gv.o frechet_gv.f

#----------------------------------

frechet_cv: $(MINEOSBIN)/frechet_cv

$(MINEOSBIN)/frechet_cv: frechet_cv.o
	-rm $(MINEOSBIN)/frechet_cv
	gfortran $(FFLAGS) -o $(MINEOSBIN)/frechet_cv frechet_cv.o

frechet_cv.o: frechet_cv.f
	gfortran $(FFLAGS) -c -o frechet_cv.o frechet_cv.f

#----------------------------------

draw_frechet_gv: $(MINEOSBIN)/draw_frechet_gv

$(MINEOSBIN)/draw_frechet_gv: draw_frechet_gv.o
	-rm $(MINEOSBIN)/draw_frechet_gv
	gfortran $(FFLAGS) -o $(MINEOSBIN)/draw_frechet_gv draw_frechet_gv.o

draw_frechet_gv.o: draw_frechet_gv.f
	gfortran $(FFLAGS) -c -o draw_frechet_gv.o draw_frechet_gv.f

#----------------------------------
clean:
	-rm *.o
