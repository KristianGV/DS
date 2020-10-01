#############################################
### Don't change anything below this line ###
#############################################


#FF=$(FC)
FF=gfortran

### Compiler options ###

# Options for linux
FOPT = -O -ffixed-line-length-none -fopenmp

### Setups for the DarkSUSY install directory ###

# Determine where to install stuff (prefix is set in configure)
prefix=/home/kriss/DS/darksusy-6.2.2
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}

LIB=$(DS_INSTALL)/lib
INC=-I./ -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include 
INC_MSSM=-I$(DS_INSTALL)/src_models/mssm/include
cfitsio=.
ADD_SCR=$(DS_INSTALL)/scr/add_libs.pl

## User libraries with user replaced functions, old just for testing
#USER_LIB=-lds_core_user -lds_$(DS_MODULE)_user
## the below is old, just for testing
#ALLLIBS= $(USER_LIB) -lds_core -lds_$(DS_MODULE)
#ifeq ($(DS_MODULE), mssm)
#  ALLLIBS += -lisajet -lHB -lFH -lisospin
#endif
## Libraries. List here the libraries and object files you want to link to
## IN THE ORDER in which you want them to be searched
#ALL_LIB = $(LIB)/libds_core_user.a $(LIB)/libds_core.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_$(DS_MODULE).a


AUX_LIB_MSSM = -lisajet -lHS -lHB -lFH -lisospin

all: dsmain_wimp dsmain_decay

#dscheckmod :
#	test -f $(DS_INSTALL)/lib/libds_${DS_MODULE}.a || { echo ERROR: Module ${DS_MODULE} does not exist, or is not compiled; exit 1;}

dscheckmod :
	@test `ls $(LIB) | grep libds_${DS_MODULE}.a` || { echo ERROR: Module ${DS_MODULE} does not exist, or is not compiled; exit 1;}

dsmain_wimp : DS_MODULE = $(shell sed -n '1p' dsmain_wimp.driver)

dsmain_wimp : dscheckmod makefile dsmain_wimp.F
dsmain_wimp : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	printf "#define MODULE_CONFIG MODULE_"$(DS_MODULE)"\n"  > module_compile.F
	printf "$(LIB)/libds_"$(DS_MODULE)"_user.a\n"$(LIB)"/libds_core_user.a\n"$(LIB)"/libds_"$(DS_MODULE)".a\n"$(LIB)"/libds_core.a" > module_link.txt
	$(ADD_SCR) libds_tmp.a module_link.txt
	$(FF) $(FOPT) $(INC) $(INC_MSSM) -L$(LIB) -o dsmain_wimp dsmain_wimp.F \
	libds_tmp.a $(shell if [ "x$(DS_MODULE)" = "xmssm" ]; then printf "%s" " $(AUX_LIB_MSSM)"; fi)
	rm -f module_compile.F
	rm -f module_link.txt
	rm -f libds_tmp.a

dsmain_decay : DS_MODULE = generic_decayingDM

dsmain_decay : dscheckmod makefile dsmain_decay.F
dsmain_decay : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	printf "#define MODULE_CONFIG MODULE_"$(DS_MODULE)"\n"  > module_compile.F
	printf "$(LIB)/libds_"$(DS_MODULE)"_user.a\n"$(LIB)"/libds_core_user.a\n"$(LIB)"/libds_"$(DS_MODULE)".a\n"$(LIB)"/libds_core.a" > module_link.txt
	$(ADD_SCR) libds_tmp.a module_link.txt
	$(FF) $(FOPT) $(INC) $(INC_MSSM) -L$(LIB) -o dsmain_decay dsmain_decay.F \
	libds_tmp.a 
	rm -f module_compile.F
	rm -f module_link.txt
	rm -f libds_tmp.a


.NOTPARALLEL: 

.PHONY: dscheckmod
