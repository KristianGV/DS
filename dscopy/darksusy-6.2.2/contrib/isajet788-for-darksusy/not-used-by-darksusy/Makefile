########################################################################
# Makefile for Isajet >=7.87 using cpp.
# Usage:
# 1. Untar isajet.tar in an empty directory.
# 2. Edit this Makefile and choose appropriate values for FC, FFLAGS,
#    and LFLAGS/LIBS. See PILOT.h for other FFLAGS options. These are
#    passed to the subdirectory Makefiles.
# 3. Type "make".
########################################################################

SHELL = /bin/csh
.IGNORE:

# Fortran
FC = gfortran
export FC

# Fortran flags -- do not want -c for link step 
# See isainc/PILOT.h for possible -D options
#FFLAGS = -cpp -O
FFLAGS = -cpp -O -DLINUX_X
#FFLAGS = -cpp -O -DMACOS_X
#FFLAGS = -cpp -O -DSUN_X
export FFLAGS

# Need both libisajet.a and libisared.a unless -DNOISATOOLS_X
#LFLAGS = -L. -lisajet
#LIBS := libisajet.a
LFLAGS = -L. -lisajet -lisared
LIBS := libisajet.a libisared.a

all: isajet.x isasusy.x isasugra.x isajet.ps

isajet.x: runjet.f $(LIBS)
	$(FC) -o isajet.x runjet.f isajet/aldata.o $(FFLAGS) -I./isainc $(LFLAGS) 

isasusy.x: ssrun.f $(LIBS)
	$(FC) -o isasusy.x ssrun.f isajet/aldata.o $(FFLAGS) -I./isainc $(LFLAGS) 

isasugra.x: sugrun.f $(LIBS)
	$(FC) -o isasugra.x sugrun.f isajet/aldata.o $(FFLAGS) -I./isainc $(LFLAGS) 

libisajet.a: 
	cd isajet && $(MAKE)
	ln -sf isajet/libisajet.a .

libisared.a: 
	cd isared && $(MAKE)
	ln -sf isared/libisared.a .

isajet/aldata.o: isajet/aldata.f
	cd isajet && $(MAKE) aldata.o

isajet.ps: isatex
	cd isatex && $(MAKE)
	ln -s isatex/isajet.ps .

clean:
	rm libisajet.a libisared.a isajet.x isasusy.x isasugra.x
	cd isajet && $(MAKE) clean
	cd isared && $(MAKE) clean
	cd isatex && $(MAKE) clean
