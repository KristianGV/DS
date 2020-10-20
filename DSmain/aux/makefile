# This folder 'aux' is for auxiliary programs that some users might
# find useful. The aim here is that both the auxiliary main programs, but also
# this makefile should prove useful when you design your own programs (and
# makefiles).
#
# Notes regarding makefiles: In DarkSUSY, we use autoconf which creates
# makefiles from templates, files with names makefile.in. ./configure then
# creates a makefile from the corresponsing makefile.in. Hence, if you change
# the makefile here, make sure to change examples/aux/makefile.in and not
# examples/aux/makefile as the latter will be overwritten everytime you
# run configure. If you copy this makefile (or makefile.in) to use as a
# starting point for your own makefiles, you can either use makefile.in
# as your starting point (if you intend to use it in an autoconf setting,
# or use makefile as your staring point if you don't use autoconf).
# During configure, when makefile is created, all variables of type
# @NAME@ in makefile.in will be replaced with their value set in the
# configure script.
#
# Compared to the examples/makefile.in, we here try to setup the
# makefile.in in a simpler way as we don't need the flexibility of changing
# particle physics modules as we do in examples/makefile.in.

#FF=$(FC)
FF=gfortran

### Compiler options ###

# Options for linux
FOPT = -O -ffixed-line-length-none -fopenmp

# DEBUG
DEBUG = -g 

#DEBUG EXTRA
DEBUG_EXTRA =  -Wextra -Wall -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan

### Setups for the DarkSUSY install directory ###

# Determine where to install stuff (prefix is set in configure)
prefix=/home/kriss/DS/darksusy-6.2.2
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}

LIB=$(DS_INSTALL)/lib
INC=-I./ -I./user_replaceables -I$(DS_INSTALL)/src/user_replaceables -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include
# These are include directories for particle physics modules (if needed)
INC_MSSM=-I$(DS_INSTALL)/src_models/mssm/include
INC_GENERIC=-I$(DS_INSTALL)/src_models/generic_wimp/include
INC_SILVEIRAZEE=-I$(DS_INSTALL)/src_models/silveira_zee/include
INC_vdSIDM=-I$(DS_INSTALL)/src_models/vdSIDM/include
cfitsio=.

# The script below creates a temporary library from different input libraries
# and object files. The purpose of this is to make sure we include
# files in the proper order (i.e. in the order they appear as arguments)
# so that default object files can be replaced by user object files.
ADD_SCR=$(DS_INSTALL)/scr/add_libs.pl

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
$(DS_INSTALL)/src/aux_xcernlib/dgadap.f $(DS_INSTALL)/src/fi/dsfi2to2oh2.f\
$(DS_INSTALL)/src_models/silveira_zee/an/dsanwx.f\
$(DS_INSTALL)/src/rd/dsrdomega.f $(DS_INSTALL)/src_models/silveira_zee/an/dssigmavpartial.f \
$(DS_INSTALL)/src/fi/dsinterpolatetable2D.f


#$(DS_INSTALL)/src/rd/dsrdthav.f 

all: scalar_FIMP abundance_mO


### We below give makefile commands to build our example programs.
### These can also serve as examples on how to write your own makefiles
### for your own programs.
### For each main program, we need to specify which particle physics module
### to use and where the include files are located for this module.
### Hence, to switch module you need to change the DS_MODULE and INC_MODULE
### declarations for that program (you might also need to add more
### libraries to include if your module relies on other packages).
### For the concept of replaceable functions to work on all platforms
### we provide a script, called via $(ADD_SCR) that adds libraries together
### and it will pick the objects in these libraries in the order of the
### arguments. Hence, in the examples below, the user replaceable objects
### (given in libds_*user.a) will replace the default ones as those libraries
### come before the default libraries in the list of libraries to add.


scalar_FIMP : DS_MODULE = silveira_zee
scalar_FIMP : INC_MODULE = -I$(DS_INSTALL)/src_models/silveira_zee/include
scalar_FIMP : scalar_FIMP.f $(DS_INSTALL)/src_models/silveira_zee/user_replaceables/dssmgammahpartial.f
scalar_FIMP : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) $(DEBUG) -o scalar_FIMP.out scalar_FIMP.f \
	$(DS_INSTALL)/src_models/silveira_zee/user_replaceables/dssmgammahpartial.f \
	libds_tmp.a
	rm -f libds_tmp.a


abundance_mO : DS_MODULE = silveira_zee
abundance_mO : INC_MODULE = -I$(DS_INSTALL)/src_models/silveira_zee/include
abundance_mO : abundance_mO.f $(DS_INSTALL)/src/user_replaceables/dsrdthav.f
abundance_mO : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) $(SRC) $(DEBUG) -o abundance_mO.out abundance_mO.f \
	$(DS_INSTALL)/src/user_replaceables/dsrdthav.f \
	libds_tmp.a
	rm -f libds_tmp.a



.NOTPARALLEL: 
