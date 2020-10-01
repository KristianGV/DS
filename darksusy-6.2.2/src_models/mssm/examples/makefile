# Makefile for src_models/mssm/examples directory
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Changed by Paolo Gondolo (2011), Joakim Edsjo (2013, 2014, 2015)
# This file is automatically created by makemf.pl.
# To aviod missing libraries this automatically created makefile.in     
# links to all libraries that could typically be needed. 
# If you use this makefile.in as a template for your own makefile
# you can of course delete the linking to not needed libraries.
#
# Define fortran compiler and options (set when ./configure is run
# in the DarkSUSY root directory).
FF=gfortran
FOPT=-O -ffixed-line-length-none -fopenmp

CC=gcc
COPT=-O

# Determine where DarkSUSY is installed
prefix=/home/kriss/darksusy-6.2.2
DS_INSTALL=${prefix}

LIB=$(DS_INSTALL)/lib
INC=-I./ -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include -I$(DS_INSTALL)/src/templates
INC_MODULE=-I$(DS_INSTALL)/src_models/mssm/include    

LIBDEPS=$(LIB)/libds_core.a $(LIB)/libds_mssm.a $(LIB)/libisajet.a
LIBS=-lds_core -lds_mssm -lds_mssm_user -lFH -lHB -lisospin -lisajet

all: dsslha2slha dstemp dstest-isasugra

dsslha2slha: dsslha2slha.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o dsslha2slha dsslha2slha.f \
	 $(LIBS)

dstemp: dstemp.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o dstemp dstemp.f \
	 $(LIBS)

dstest-isasugra: dstest-isasugra.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o dstest-isasugra dstest-isasugra.f \
	 $(LIBS)

