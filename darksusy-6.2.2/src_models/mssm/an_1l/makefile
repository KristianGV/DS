# Makefile for src_models/mssm/an_1l directory
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

INC_DEP = dsidtag.h dsmpconst.h dsmssm.h dsparticles.h dssm.h

vpath %.h $(INC)

SRC = dsanggim.f dsanggimpar.f dsanggre.f dsanggrepar.f dsanglglim.f \
dsanglglre.f dsanzg.f dsanzgpar.f dsfl1c1.f dsfl1c2.f dsfl2c1.f \
dsfl2c2.f dsfl3c1.f dsfl3c2.f dsfl4c1.f dsfl4c2.f dsi_12.f dsi_13.f \
dsi_14.f dsi_22.f dsi_23.f dsi_24.f dsi_32.f dsi_33.f dsi_34.f dsi_41.f \
dsi_42.f dsilp2.f dsj_1.f dsj_2.f dsj_3.f dslp2.f dspi1.f dspiw2.f \
dspiw2i.f dspiw3.f dspiw3i.f dsrepfbox.f dsrepgh.f dsrepw.f dsslc1.f \
dsslc2.f dssubka.f dssubkb.f dssubkc.f dsti_214.f dsti_224.f dsti_23.f \
dsti_33.f dsti_5.f

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
