# Makefile for src/se_mod directory
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

INC_DEP = dsio.h dsmpconst.h dssecom.h dssem_earth.h dssem_sun.h

vpath %.h $(INC)

SRC = dssem_earthdens.f dssem_earthdenscomp.f dssem_earthmass.f \
dssem_earthmassint.f dssem_earthne.f dssem_earthpot.f \
dssem_earthpotint.f dssem_earthvesc.f dssem_edfunc.f dssem_epfunc.f \
dssem_spfunc.f dssem_suncdens.f dssem_suncdensint.f dssem_suncdfunc.f \
dssem_sundens.f dssem_sundenscomp.f dssem_sunmass.f dssem_sunmfrac.f \
dssem_sunne.f dssem_sunne2x.f dssem_sunpot.f dssem_sunpotint.f \
dssem_sunread.f dssem_sunset.f dssem_sunvesc.f dssem_sunx2z.f \
dssem_sunz2x.f

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
