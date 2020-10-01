# Makefile for src_models/mssm/ac_bsg directory
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Changed by Paolo Gondolo (2011), Joakim Edsjo (2013)
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
INC=../include ../../include ../../../src/include ../../../contrib/include
DINC=-I../include -I../../include -I../../../src/include -I../../../contrib/include
DOBJ=../../../tmp/build-mssm

INC_DEP = dsmpconst.h dsmssm.h dsparticles.h dssm.h

vpath %.h $(INC)

SRC = dsbsg2007full.f dsbsgalpha3.f dsbsgalpha3int.f dsbsgammafull.f \
dsbsgat0.f dsbsgat1.f dsbsgbofe.f dsbsgc41h2.f dsbsgc41susy.f \
dsbsgc70h2.f dsbsgc70susy.f dsbsgc71chisusy.f dsbsgc71h2.f \
dsbsgc71phi1susy.f dsbsgc71phi2susy.f dsbsgc71wsusy.f dsbsgc80h2.f \
dsbsgc80susy.f dsbsgc81chisusy.f dsbsgc81h2.f dsbsgc81phi1susy.f \
dsbsgc81phi2susy.f dsbsgc81wsusy.f dsbsgckm.f dsbsgd1td.f dsbsgd2d.f \
dsbsgd2td.f dsbsgd7chi1.f dsbsgd7chi2.f dsbsgd7h.f dsbsgd8chi1.f \
dsbsgd8chi2.f dsbsgd8h.f dsbsgeb.f dsbsgechi.f dsbsgeh.f dsbsget0.f \
dsbsgf71.f dsbsgf72.f dsbsgf73.f dsbsgf81.f dsbsgf82.f dsbsgf83.f \
dsbsgft0.f dsbsgft1.f dsbsgfxy.f dsbsgg.f dsbsgg7chi2.f dsbsgg7chij1.f \
dsbsgg7h.f dsbsgg7w.f dsbsgg8chi2.f dsbsgg8chij1.f dsbsgg8h.f dsbsgg8w.f \
dsbsggxy.f dsbsgh1x.f dsbsgh2xy.f dsbsgh3x.f dsbsghd.f dsbsghtd.f \
dsbsgkc.f dsbsgkt.f dsbsgmtmuw.f dsbsgmtmuwint.f dsbsgpre2007full.f \
dsbsgri.f dsbsgud.f dsbsgutd.f dsbsgwud.f dsbsgwxy.f dsbsgyt.f \
dsphi22a.f dsphi22b.f dsphi27a.f dsphi27b.f

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
