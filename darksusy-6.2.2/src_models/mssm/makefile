# Makefile for the DarkSUSY library
# Technical questions about this makefile should be directed to
# Joakim Edsjo, edsjo@fysik.su.se
# Date: August 29, 2000.

# You shouldn't have to change anything in this file.
# You should only change in the makefile in the DarkSUSY root directory.

###############################################
#### DO NOT CHANGE ANYTHING BELOW THIS LINE ###
###############################################

# Here we set where the library goes and where include files are.

LIB=../../lib
#TSPACKDIR=/home/kriss/darksusy-6.2.2/contrib/TSPACK
BUILD=../../tmp/build-mssm
BUILDCOMMON=../../tmp/build-common

vpath %.a $(LIB)

DIRS = ac ac_bsg an an_1l an_ib an_ib2 an_ib3 an_sf an_stu an_yield \
an_yield_casc cr dd examples_aux ge ge_cmssm ge_slha ini kd rd se_yield \
se_yield_casc xfeynhiggs xhiggsbounds xsuperiso ../common/aux \
../common/sm

USER_DIRS = user_replaceables

DOBJ = $(DIRS:%=%_obj)

USER_DOBJ = $(USER_DIRS:%=%_obj)

# Trick to call main makefile in case make is called in src directory.
all-local : make-lib

make-lib :
	cd ../..; make ds_mssm

all : libds_mssm_user.a libds_mssm.a 

libds_mssm.a : $(DOBJ)
	ar rS $(LIB)/libds_mssm.a $(BUILDCOMMON)/*.o
	ar rS $(LIB)/libds_mssm.a $(BUILD)/*.o
	ranlib $(LIB)/libds_mssm.a

libds_mssm_user.a: $(USER_DOBJ)

$(DOBJ) : 
	cd $(@:_obj=) ; make all

$(USER_DOBJ) : 
	cd $(@:_obj=) ; make all

.PHONY: clean distclean $(DOBJ) $(USER_DOBJ)
distclean : clean
	rm $(LIB)/libds_mssm.a
	rm $(LIB)/libds_mssm_user.a

clean : 
	rm -f */*.o 

