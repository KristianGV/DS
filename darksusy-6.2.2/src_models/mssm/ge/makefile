# Makefile for src_models/mssm/ge directory
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

INC_DEP = dsidtag.h dsio.h dsisasugra.h dsmpconst.h dsmssm.h \
dsparticles.h dssm.h

vpath %.h $(INC)

SRC = dsb0loop.f dschasct.f dsdmspin.f dsg0loop.f dsg4set.f dsg4set12.f \
dsg4set1234.f dsg4set13.f dsg4set23.f dsg4set34.f dsg4setc.f \
dsg4setc12.f dsg4setc1234.f dsg4setc13.f dsg4setc23.f dsg4setc34.f \
dshgfu.f dshigsct.f dshigwid.f dshigwid1.f dshlf2.f dshlf3.f dsmass.f \
dsmwimp.f dsneusct.f dsorder_flavour.f dspole.f dsqindx.f dsralph3.f \
dsrghm.f dsrmq.f dssfesct.f dsspectrum.f dsspwid.f dssuconst_higgs.f \
dssuconst_yukawa.f dssuconst_yukawa_running.f dsvertx.f dsvertx1.f \
dsvertx3.f dswhwarn.f dswspectrum.f dswunph.f dswvertx.f dswwidth.f \
g4p.f

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
