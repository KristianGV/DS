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
BUILD=../../tmp/build-generic_wimp
BUILDCOMMON=../../tmp/build-common

vpath %.a $(LIB)

# Here we list all directories which contains files to be made
DIRS = an cr dd ge ini kd rd se_yield ../common/aux ../common/sm

USER_DIRS = user_replaceables

DOBJ = $(DIRS:%=%_obj)

USER_DOBJ = $(USER_DIRS:%=%_obj)

# Trick to call main makefile in case make is called in src directory.
all-local : make-lib

make-lib :
	cd ../..; make ds_generic_wimp

all : libds_generic_wimp_user.a libds_generic_wimp.a 

libds_generic_wimp.a : $(DOBJ)
	ar rS $(LIB)/libds_generic_wimp.a $(BUILDCOMMON)/*.o
	ar rS $(LIB)/libds_generic_wimp.a $(BUILD)/*.o
	ranlib $(LIB)/libds_generic_wimp.a

libds_generic_wimp_user.a: $(USER_DOBJ)

$(DOBJ) : 
	cd $(@:_obj=) ; make all

$(USER_DOBJ) : 
	cd $(@:_obj=) ; make all

.PHONY: clean distclean $(DOBJ) $(USER_DOBJ)
distclean : clean
	rm $(LIB)/libds_generic_wimp.a
	rm $(LIB)/libds_generic_wimp_user.a

clean : 
	rm -f */*.o 

