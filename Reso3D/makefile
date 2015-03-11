# makefile for the azhydro0p2 programs
#   P.K. September 2004
#   E.F  September 2005
#
CC = gcc
CFLAGS =  $(GRAPHINC) -O3
NAGLIB = /psi/nagfl-18/lib/libnag.a /psi/NAGWare_f95-4.0/lib/libf95.a
#NAGLIB=/psi/nagfl-17d9/lib/libf90.a
#
# targets
# run make -f makefile azhydro0p2 to create the hydro program.
#--------------- hydro evolution ------------------
azhydro0p2.o:		azhydro0p2.f
			gfortran -ff2c -O -c azhydro0p2.f

azinitials0p2.o:	azinitials0p2.f
			gfortran -ff2c -O -c azinitials0p2.f

azhyspeca0p2.o:		azhyspeca0p2.f
			gfortran -ff2c -O -c azhyspeca0p2.f

azhycom0p2.o:		azhycom0p2.f
			gfortran -ff2c -O -c azhycom0p2.f

azshastla0p2.o:		azshastla0p2.f
			gfortran -ff2c -O -c azshastla0p2.f

azhydro0p2:	azhydro0p2.o azinitials0p2.o azhyspeca0p2.o azhycom0p2.o  azshastla0p2.o
		gfortran -ff2c -O azhydro0p2.o azinitials0p2.o azhyspeca0p2.o azhycom0p2.o azshastla0p2.o  -o azhydro0p2 

#--------------- calculation of spectra-------------------------------

azspectra0p2.o:	azspectra0p2.f
		gfortran -ff2c -O -c azspectra0p2.f

azbesselK0p2.o:	azbesselK0p2.f
		gfortran -ff2c -O -c azbesselK0p2.f

azspectra0p2:	azspectra0p2.o azhycom0p2.o azbesselK0p2.o azinitials0p2.o
		gfortran -ff2c -O azspectra0p2.o azhycom0p2.o azbesselK0p2.o azinitials0p2.o -o azspectra0p2


#************** calculation of resonance decays ****************

functions.o:	functions.c functions.h reso.h
		$(CC) $(CFLAGS)  -c functions.c

decay.o:	decay.c decay.h
		$(CC) $(CFLAGS) -c decay.c  

int.o:		int.c int.h 
		$(CC) $(CFLAGS) -c int.c  

tools.o:	tools.c tools.h 
		$(CC) $(CFLAGS) -c tools.c  

reso.o:		reso.c reso.h functions.h 
		$(CC) $(CFLAGS) -c reso.c  

reso:		reso.o functions.o int.o tools.o decay.o reso.h functions.h
		$(CC) $(CFLAGS) -o reso reso.o functions.o decay.o int.o tools.o -lm 

#-lgsl -lgslcblas

