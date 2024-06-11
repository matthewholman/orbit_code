# Make file for orbit fitting programs.  Change flags below as required.

include Makefile.defs

CC = gcc
CFLAGS = -std=c99 -O -Wall -fPIC
LIBS = -lm -lrebound -lassist

PROGS = fit_radec
ALL:  $(PROGS)

clean:
	rm *.o $(PROGS)

OBJS  = fit_radec.o

NR   = nrutil.o ludcmp.o lubksb.o gaussj.o mrqmin_orbit.o mrqcof_orbit.o \
	covsrt.o gasdev.o ran1.o
#ORBSUBS = orbfit1.o transforms.o ephem_earth.o
ORBSUBS = orbfit1.o transforms.o


//ephem_earth.o: ephem_read.h
//orbfit1.o ephem_earth.o: ephem_types.h
$(OBJS) $(ORBSUBS) orbfit2.o mrqmin_orbit.o mrqcof_orbit.o :  orbfit.h Makefile

fit_radec: fit_radec.o orbfit2.o $(ORBSUBS) $(NR)
	$(CC) fit_radec.o orbfit2.o $(ORBSUBS) $(NR) $(LIBS) -L. -o fit_radec

cloud: cloud.o orbfit2.o $(ORBSUBS) $(NR)
	$(CC) cloud.o orbfit2.o $(ORBSUBS) $(NR) $(LIBS) -L. -o cloud

fit_radec_prelim: fit_radec_prelim.o orbfit2.o $(ORBSUBS) $(NR)
	$(CC) fit_radec_prelim.o orbfit2.o $(ORBSUBS) $(NR) $(LIBS) -o fit_radec_prelim

library: orbfit2.o $(ORBSUBS) $(NR)
	$(CC) -shared -fPIC orbfit2.o $(ORBSUBS) $(NR) $(LIBS) -o liborbfit.so
