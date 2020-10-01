# Makefile for src/dd directory
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Changed by Paolo Gondolo (2011), Joakim Edsjo (2013, 2014, 2015, 2016, 2017)
# This file is automatically created by makemf.pl.

# Define fortran compiler and options (set when ./configure is run
# in the DarkSUSY root directory).
FF=gfortran
FOPT=-O -ffixed-line-length-none -fopenmp

FC=$(FF)
FFLAGS=$(FOPT) -c $(DINC)

CC=gcc
COPT=-O
CCFLAGS=$(COPT) -c $(DINC)

# Dependencies and libraries
INC=../include ../../contrib/include
DINC=-I../include -I../../contrib/include
DOBJ=../../tmp/build-src

INC_DEP = dsddcom.h dshmcom.h dsio.h dsmpconst.h dsnuclides.h

vpath %.h $(INC)

SRC = dsddDMCRcountrate.f dsddDMCRdgammadt.f dsddDMCRflux.f \
dsddDMCRquenching.f dsddDMCRquenching_set.f dsddDMCRsigCR.f \
dsddDMCRsigtarget.f dsddTDMattenuation.f dsddTrmax.f dsdddrde.f \
dsddeta.f dsddfffb.f dsddfffermi.f dsddffgauss.f dsddffgould.f \
dsddffh41.f dsddffism.f dsddffls.f dsddffoddg.f dsddffsd.f dsddffsi.f \
dsddffsimsd.f dsddffsog.f dsddffspsm.f dsddg2sigma.f dsddhelp.f \
dsddinit.f dsddlfreesimp.f dsddlim.f dsddlimits.f dsddpp.f dsddreadQF.f \
dsddset.f dsddsigmanucleon.f dsddvearth.f dsnucldindx.f dsz_of.f

OBJ1 = $(patsubst %.f,$(DOBJ)/%.o,$(SRC))

OBJ2 = $(patsubst %.F,$(DOBJ)/%.o,$(OBJ1))

OBJ3 = $(patsubst %.f90,$(DOBJ)/%.o,$(OBJ2))

OBJ4 = $(patsubst %.F90,$(DOBJ)/%.o,$(OBJ3))

OBJ = $(patsubst %.c,$(DOBJ)/%.o,$(OBJ4))

all : $(OBJ)  mods

$(DOBJ)/%.o : %.F $(INC_DEP) 
	$(FC) $(FFLAGS) $< -o $@

$(DOBJ)/%.o : %.f $(INC_DEP) 
	$(FC) $(FFLAGS) $< -o $@

$(DOBJ)/%.o : %.F90 $(INC_DEP) 
	$(FC) $(FFLAGS) $< -o $@

$(DOBJ)/%.o : %.f90 $(INC_DEP) 
	$(FC) $(FFLAGS) $< -o $@

$(DOBJ)/%.o : %.c $(INC_DEP) 
	$(CC) $(CCFLAGS) $< -o $@


mods : 
	if [ -f *.mod ]; then\
		cp -p *.mod ../include/ ;\
	fi
