# Makefile for src_models/mssm/an_ib2 directory
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

INC_DEP = dsanyieldcom.h dsanyieldmodelcom.h dsib2com.h dsio.h \
dsmpconst.h dsmssm.h dsparticles.h dssm.h

vpath %.h $(INC)

SRC = dsib2BRtree.f dsib2Msqaux.f dsib2chinit.f dsib2convint.f \
dsib2dnde.f dsib2dsde_aux.f dsib2dsde_auxfsr.f dsib2ffA0amps.f \
dsib2ffA0sFSR.f dsib2ffA0sISR.f dsib2ffA0sVIB.f dsib2ffA0tuFSR.f \
dsib2ffA0tuISR.f dsib2ffA0tuVIB.f dsib2ffH0amps.f dsib2ffH0sFSR.f \
dsib2ffH0sISR.f dsib2ffH0sVIB.f dsib2ffH0tuFSR.f dsib2ffH0tuISR.f \
dsib2ffH0tuVIB.f dsib2ffHPamps.f dsib2ffHPsFSRH3.f dsib2ffHPsFSRZ.f \
dsib2ffHPsISRH.f dsib2ffHPsISRW.f dsib2ffHPsVIB.f dsib2ffHPtuFSR.f \
dsib2ffHPtuISR1.f dsib2ffHPtuISR2.f dsib2ffHPtuVIB.f dsib2ffWamps.f \
dsib2ffWsFSR.f dsib2ffWsISR.f dsib2ffWsVIB.f dsib2ffWtuFSR.f \
dsib2ffWtuISR1.f dsib2ffWtuISR2.f dsib2ffWtuVIB.f dsib2ffZamps.f \
dsib2ffZsFSR.f dsib2ffZsISR.f dsib2ffZsVIB.f dsib2ffZtuFSR.f \
dsib2ffZtuISR.f dsib2ffZtuVIB.f dsib2getsfermion.f dsib2intres.f \
dsib2kinematics.f dsib2set.f dsib2sigmav.f dsib2yield.f dsib2yieldone.f

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
