#FF=$(FC)
FF=gfortran

### Compiler options ###

# Options for linux
FOPT = -O -ffixed-line-length-none -fopenmp

### Setups for the DarkSUSY install directory ###

# Determine where to install stuff (prefix is set in configure)
prefix=/home/kriss/darksusy-6.2.2
# DS_INSTALL is where the library and data files will be installed
DS_INSTALL=${prefix}

LIB=$(DS_INSTALL)/lib
INC=-I./ -I$(DS_INSTALL)/src/include -I$(DS_INSTALL)/src_models/include -I$(DS_INSTALL)/contrib/include
INC_MSSM=-I$(DS_INSTALL)/src_models/mssm/include
cfitsio=.
ADD_SCR=$(DS_INSTALL)/scr/add_libs.pl


# User libraries with user replaced functions (this only works on OS X, and is
# NOT the default building method anymore)
#USER_LIB=-lds_core_user -lds_$(DS_MODULE)_user
# Libraries. List here the libraries and object files you want to link to
# IN THE ORDER in which you want them to be searched
#ALL_LIB = $(LIB)/libds_core_user.a $(LIB)/libds_core.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_$(DS_MODULE).a


all: dstest_mssm dstest_silveira_zee dstest_genWIMP

dstest_mssm : DS_MODULE = mssm
dstest_mssm : dstest_mssm.f $(LIB)/libisajet.a $(LIB)/libHS.a $(LIB)/libHB.a $(LIB)/libFH.a
dstest_mssm : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	printf "#define MODULE_CONFIG MODULE_"$(DS_MODULE)"\n"  > module_compile.F
	printf "$(LIB)/libds_"$(DS_MODULE)"_user.a\n"$(LIB)"/lib/ds_core_user.a\n"$(LIB)"/libds_"$(DS_MODULE)".a\n"$(LIB)"/libds_core.a" > module_link.txt
	$(ADD_SCR) libds_tmp.a module_link.txt
	$(FF) $(FOPT) $(INC) $(INC_MSSM) -L$(LIB) -o dstest_mssm dstest_mssm.f \
	libds_tmp.a -lHS -lHB -lFH -lisospin -lisajet
	rm -f module_compile.F
	rm -f module_link.txt
	rm -f libds_tmp.a

dstest_shared : DS_MODULE = mssm
dstest_shared : dstest_mssm.f $(LIB)/libds_core_mssm.so
	printf "#define MODULE_CONFIG MODULE_"$(DS_MODULE)"\n"  > module_compile.F
	$(FF) $(FOPT) $(INC) $(INC_MSSM) -L$(LIB) -o dstest_shared dstest_mssm.f \
	-lds_core_mssm
	rm -f module_compile.F

dstest_silveira_zee : DS_MODULE = silveira_zee
dstest_silveira_zee : INC_MODULE = -I$(DS_INSTALL)/src_models/silveira_zee/include
dstest_silveira_zee : dstest_silveira_zee.f
dstest_silveira_zee : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o dstest_silveira_zee dstest_silveira_zee.f \
	libds_tmp.a 
	rm -f libds_tmp.a

dstest_genWIMP : DS_MODULE = generic_wimp
dstest_genWIMP : INC_MODULE = -I$(DS_INSTALL)/src_models/silveira_zee/include
dstest_genWIMP : dstest_genWIMP.f
dstest_genWIMP : $(LIB)/libds_core.a $(LIB)/libds_core_user.a
	$(ADD_SCR) libds_tmp.a $(LIB)/libds_$(DS_MODULE)_user.a $(LIB)/libds_core_user.a $(LIB)/libds_$(DS_MODULE).a $(LIB)/libds_core.a
	$(FF) $(FOPT) $(INC) $(INC_MODULE) -L$(LIB) -o dstest_genWIMP dstest_genWIMP.f \
	libds_tmp.a 
	rm -f libds_tmp.a

.NOTPARALLEL: 
