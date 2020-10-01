# Makefile for src/fi directory
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
INC=-I./ -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include


SRC = $(DS_INSTALL)/src/aux/dsbessei1.f \
$(DS_INSTALL)/src/aux/dsbessek1.f $(DS_INSTALL)/src/aux/dsisnan.f \
$(DS_INSTALL)/src/aux/dsf_int.f $(DS_INSTALL)/src/aux/dsf_int2.f \
$(DS_INSTALL)/src/aux/dswrite.f $(DS_INSTALL)/src/aux/dsi_trim.f \
$(DS_INSTALL)/src/rd/dsrddof.f $(DS_INSTALL)/src/kd/dskdgeff.f \
$(DS_INSTALL)/src/rd/dsrdset.f $(DS_INSTALL)/src/rd/dsrdreaddof.f \
$(DS_INSTALL)/src/aux/dsdatafile.f $(DS_INSTALL)/src/aux/dsidtag.f \
$(DS_INSTALL)/src/fi/dsfidecint.f $(DS_INSTALL)/src/fi/dsfistat.f  \
$(DS_INSTALL)/src/fi/dsfik1.f $(DS_INSTALL)/src/fi/dsfipcm.f \
$(DS_INSTALL)/src/fi/dsfidecab.f $(DS_INSTALL)/src/fi/dsfidecoh2.f\
$(DS_INSTALL)/src/fi/dsfi2to2int.f $(DS_INSTALL)/src/fi/dsfi2to2rhs.f\
$(DS_INSTALL)/src/fi/dsfisigma.f $(DS_INSTALL)/src/fi/dsfi2to2ab.f\
$(DS_INSTALL)/src/aux_xcernlib/dgadap.f $(DS_INSTALL)/src/fi/dsfi2to2oh2.f

all: test dydt 2to2rhs test2to2 generic_FIMP

test: test.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) -o test test.f \
	 $(LIBS)

test2to2: test2to2.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) -o test2to2 test2to2.f \
	 $(LIBS)

dydt: dydt.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) -o dydt dydt.f \
	 $(LIBS)

2to2rhs: 2to2rhs.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) -o 2to2rhs 2to2rhs.f \
	 $(LIBS)

generic_FIMP: generic_FIMP.f $(LIBDEPS)
	 $(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) -o generic_FIMP generic_FIMP.f \
	 $(LIBS)

