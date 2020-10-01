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
BUILD=../../tmp/build-empty
BUILDCOMMON=../../tmp/build-common

vpath %.a $(LIB)

# Here we list all directories which contains files to be made
DIRS = ac an cr dd ge ini kd rd se_yield si ../common/aux ../common/sm

USER_DIRS = user_replaceables

DOBJ = $(DIRS:%=%_obj)

USER_DOBJ = $(USER_DIRS:%=%_obj)

# Trick to call main makefile in case make is called in src directory.
all-local : make-lib

make-lib :
	cd ../..; make ds_empty

all : libds_empty_user.a libds_empty.a

libds_empty.a : $(DOBJ)
	ar rS $(LIB)/libds_empty.a $(BUILDCOMMON)/*.o
	ar rS $(LIB)/libds_empty.a $(BUILD)/*.o
	ranlib $(LIB)/libds_empty.a

libds_empty_user.a: $(USER_DOBJ)

$(DOBJ) : 
	cd $(@:_obj=) ; make all

$(USER_DOBJ) : 
	cd $(@:_obj=) ; make all

.PHONY: clean distclean $(DOBJ) $(USER_DOBJ)
distclean : clean
	rm $(LIB)/libds_empty.a
	rm $(LIB)/libds_empty_user.a

clean : 
	rm -f */*.o 


