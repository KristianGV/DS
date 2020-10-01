# Makefile for src_models/mssm/an_stu directory
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

INC_DEP = dsandiacom.h dsmssm.h dsparticles.h dssm.h

vpath %.h $(INC)

SRC = dsansffsff.f dsansffsss.f dsansffssv.f dsansffsvs.f dsansffsvv.f \
dsansffvff.f dsansffvss.f dsansffvsv.f dsansffvvs.f dsansffvvv.f \
dsantfffss.f dsantfffssex.f dsantfffssin.f dsantfffsv.f dsantfffsvin.f \
dsantfffvs.f dsantfffvsex.f dsantfffvv.f dsantfffvvex.f dsantfffvvin.f \
dsantffsff.f dsanufffss.f dsanufffssin.f dsanufffsv.f dsanufffsvin.f \
dsanufffvs.f dsanufffvv.f dsanufffvvex.f dsanufffvvin.f dsanuffsff.f

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
