#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(DarkSUSY,, edsjo@fysik.su.se)

AC_PREFIX_DEFAULT(${PWD})
#AC_CONFIG_SRCDIR([contrib/cfitsio/speed.c])
AC_CONFIG_SRCDIR([src/dd/dsdddrde.f])
#AC_CONFIG_HEADER([config.h])

# Checks for programs.
#AC_PROG_AWK
#AC_PROG_YACC
#AC_PROG_CXX
#AC_PROG_CC
#AC_PROG_LEX
AC_PROG_LN_S
#AC_PROG_MAKE_SET
AC_PROG_RANLIB

# F77 compiler
FCSET=false
if test "${FCFLAGS+set}" = set; then
   FCSET=true
fi

AC_PROG_FC([ifort gfortran])
if test "${FCSET}" = "false" ; then
  if test "${FC}" = "ifort"; then 
    echo "Fortran compiler is $FC"
    FCFLAGS=${FCFLAGS+-O -extend_source -openmp}
  fi 
  if test "${FC}" = "gfortran"; then
    echo "Fortran compiler is $FC"
    FCFLAGS=${FCFLAGS+-O -ffixed-line-length-none -fopenmp}
  fi
fi

# C compiler
CSET=false
if test "${CFLAGS+set}" = set; then
   CSET=true
fi

AC_PROG_CC([icc gcc])
AC_PROG_CXX([icc g++])
if test "${CSET}" = "false" ; then
  if test "${CC}" = "icc"; then 
    echo "C compiler is $CC"
    CFLAGS=${CFLAGS+-O}
  fi 
  if test "${CC}" = "gcc"; then
    echo "C compiler is $CC"
    CFLAGS=${CFLAGS+-O}
  fi
fi

AC_PROG_F77([ifort gfortran])
F77=${FC}
FFLAGS=${FCFLAGS}
#FFLAGS=${FFLAGS--funroll-loops -fno-emulate-complex -ffast-math -O}
# Set default g77 options
AC_ARG_VAR(FOPT,[Fortran compiler flags for DarkSUSY])
# If not set, set to default for g77
#FOPT=${FOPT--O -ffixed-line-length-none}
#FOPT=${FOPT--funroll-loops -fno-emulate-complex -ffast-math -O -ffixed-line-length-none}
#FOPT="-funroll-loops -fno-emulate-complex"
FOPT=${FCFLAGS}

# root directory
AC_ARG_VAR(DSINSTALL,[Full path to install directory])
if test "${prefix}" = "NONE"; then
DSINSTALL="${PWD}"
else
DSINSTALL="${prefix}"
fi

# ds version
AC_ARG_VAR(DSVERSION,[String containing the DarkSUSY version])
#DSVERSION=`echo $DSINSTALL | sed s:[^/]*/::g` # does not work
#DSVERSION=$DSINSTALL
DSVERSION="${PWD}"

# svn revision
AC_ARG_VAR(DSREV,[String containing DarkSUSY svn revision])
AC_CHECK_PROGS([SVN],[svn],[:])
if test "$SVN" = :; then
DSREV=""
else
DSREV=" rev ("`svn info | grep Revision | sed "s/Revision: //"`")"
fi

# Make dsdirver.h from DSVERSION and DSINSTALL
echo 'Creating include/dsdirver.h'
scr/config2.pl ${DSVERSION} ${DSINSTALL}

#####
##### NOTE FOR DarkSUSY Developers
##### When contributed packages are updated, add the correct path to them
##### below and then run 'autoconf' in the DarkSUSY
##### root to update the configure script. Then run ./configure and make
##### as usual to create makefiles and make DarkSUSY.
##### As of 2010-03-05, we no longer use symbolic links for contributed
##### packages

# ISAJET
AC_ARG_VAR(ISAJETDIR,[Directory where ISAJET is])
# Set ISAJETDIR if not set
ISAJETDIR=${ISAJETDIR-${PWD}/contrib/isajet788-for-darksusy}

# HIGGSBOUNDSDIR
AC_ARG_VAR(HIGGSBOUNDSDIR,[Directory where HiggsBounds is])
# Set HIGGSBOUNDSDIR if not set
HIGGSBOUNDSDIR=${HIGGSBOUNDSDIR-${PWD}/contrib/HiggsBounds-4.3.1}
#HIGGSBOUNDSDIR=${HIGGSBOUNDSDIR-${PWD}/contrib/HiggsBounds-3.8.1/HiggsBounds}
#HIGGSBOUNDSDIR=${HIGGSBOUNDSDIR-${PWD}/contrib/HiggsBounds-3.7.0/HiggsBounds-f90}
#HIGGSBOUNDSDIR=${HIGGSBOUNDSDIR-${PWD}/contrib/HiggsBounds-3.5.0beta-for-darksusy/HiggsBounds-f90}
AC_ARG_VAR(HBFHINCLUDE,[FeynHiggs include directory for HiggsBounds])
HBFHINCLUDE=${HBFHINCLUDE-${PWD}/include}
AC_ARG_VAR(HBFHLIBS,[FeynHiggs lib directory for HiggsBounds])
HBFHLIBS=${HBFHLIBS-${PWD}/lib}
AC_ARG_VAR(HBCPSUPERHLIBS,[CPSuperH lib directory for HiggsBounds])
HBCPSUPERHLIBS=${HBCPSUPERHLIBS-${PWD}/lib}
# Run configure in HIGGSBOUNDSDIR
echo "Running configure in ${HIGGSBOUNDSDIR}"
(cd ${HIGGSBOUNDSDIR} ; 
echo "./configure_darksusy" ;
./configure_darksusy)

# HiggsSignals
AC_ARG_VAR(HIGGSSIGNALSDIR,[Directory where HiggsSignals is])
HIGGSSIGNALSDIR=${HIGGSSIGNALSDIR-${PWD}/contrib/HiggsSignals-1.4.0}

# Run configure in HIGGSSIGNALSDIR
echo "Running configure in ${HIGGSSIGNALSDIR}"
(cd ${HIGGSSIGNALSDIR} ; 
echo "./configure_darksusy" ;
./configure_darksusy)

# SLHADIR
AC_ARG_VAR(SLHADIR,[Directory where SLHALIB is])
# Set SLHADIR if not set
SLHADIR=${SLHADIR-${PWD}/contrib/SLHALib-2.2}

# FeynHiggs
AC_ARG_VAR(FEYNHIGGSDIR,[Directory where FeynHiggs is])
# Set FeynHiggsdir if not set
FEYNHIGGSDIR=${FEYNHIGGSDIR-${PWD}/contrib/FeynHiggs-2.13.0}

# SuperIso
AC_ARG_VAR(SUPERISODIR,[Directory where SuperISO is])
# Set SuperIsoDir if not set
SUPERISODIR=${SUPERISODIR-${PWD}/contrib/superiso_v3.4}

# Run configure in FeynHiggsdir
echo "Running configure in ${FEYNHIGGSDIR}"
(cd ${FEYNHIGGSDIR} ; 
#if test "${F77}" = "g77" ; then
#   COPT="--prefix=ds FC=${F77}"
#else
#   COPT="--prefix=ds FFLAGS=default FC=${F77}"
#fi
# there is a bounds-check problem when linking FH to DS...........
echo "./configure --prefix=ds FC=${FC}" ;
./configure FC=${FC} --prefix=ds)
#echo "./configure --prefix=ds FC=${FC} FFLAGS=""${FCFLAGS}""" ;
#./configure FC=${FC} --prefix=ds FFLAGS="${FCFLAGS}")

# HealPix
AC_ARG_VAR(HEALPIXDIR,[Directory where HealPix is])
# Set HEALPIXDIR if not set
HEALPIXDIR=${HEALPIXDIR-${PWD}/contrib/Healpix_3.40}

# Run configure in HealPix
# Don't do this here, instead create configure_darkusy that is run at make time

# TSPACK
AC_ARG_VAR(TSPACKDIR,[Directory where TSPACK is])
# Set TSPACKDIR if not set
TSPACKDIR=${TSPACKDIR-${PWD}/contrib/TSPACK}


# Run configure in cfitsio
echo "Running configure in contrib/cfitsio"
(cd contrib/cfitsio; ./configure)

SLHADIR=${SLHADIR-${PWD}/contrib/slha}
# Run configure in SLHA
echo "Running configure in ${SLHADIR}"
(cd ${SLHADIR} ; 
echo "./configure --prefix=ds FC=${F77}" ;
unset FFLAGS
./configure FC=${F77} --prefix=ds)

AC_ARG_VAR(CFITSIO,[Location of cfitsio library])

# Checks for libraries (external of DarkSUSY).
AC_ARG_WITH([cfitsio],[  --with-cfitsio=CFITSIO  location of cfitsio library],[CFITSIO=$withval],[CFITSIO=.])
AC_CHECK_LIB([cfitsio], [ffopen],[rm -f .no_cfitsio],[touch .no_cfitsio],[-L$CFITSIO -lm])
# AC_CHECK_LIB([darksusy], [dsinit])
AC_CHECK_LIB([g2c], [erfc])
#AC_CHECK_LIB([ink], [main])
AC_CHECK_LIB([isajet], [sugra])
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
#AC_FUNC_ALLOCA
#AC_HEADER_STDC
#AC_CHECK_HEADERS([arpa/inet.h fcntl.h limits.h malloc.h memory.h netdb.h netinet/in.h stddef.h stdlib.h string.h sys/socket.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_CONST
#AC_TYPE_MODE_T
#AC_TYPE_SIZE_T
#AC_STRUCT_TM

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_FSEEKO
#AC_FUNC_LSTAT
#AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
#AC_FUNC_MALLOC
#AC_FUNC_MEMCMP
#AC_FUNC_REALLOC
#AC_TYPE_SIGNAL
#AC_FUNC_STRFTIME
#AC_FUNC_STRTOD
#AC_CHECK_FUNCS([alarm atexit floor ftruncate getcwd gethostbyname memchr memset pow socket sqrt strcasecmp strchr strcspn strncasecmp strrchr strspn strstr strtol])

# do not put contrib/feynhiggs/makefile here - it has its own configure
AC_CONFIG_FILES([${ISAJETDIR}/makefile
                 ${HIGGSBOUNDSDIR}/makefile_darksusy
                 ${HIGGSSIGNALSDIR}/makefile_darksusy
		 ${SUPERISODIR}/makefile
		 ${TSPACKDIR}/makefile
		 ${HEALPIXDIR}/configure_darksusy
                 data/makefile
                 makefile
                 examples/makefile
		 examples/aux/makefile
                 examples/test/makefile
                 src/an_yield/makefile
                 src/aux/makefile
                 src/aux_xcernlib/makefile
                 src/aux_xdiag/makefile
                 src/aux_xquadpack/makefile
                 src/cr_aux/makefile
                 src/cr_axi/makefile
                 src/cr_gamma/makefile
                 src/cr_nu/makefile
                 src/cr_ps/makefile
                 src/dd/makefile
                 src/dmd_astro/makefile
                 src/dmd_aux/makefile
                 src/dmd_mod/makefile
                 src/dmd_vel/makefile
                 src/ini/makefile
                 src/kd/makefile
                 src/rd/makefile
                 src/se_aux/makefile
                 src/se_mod/makefile
                 src/se_nu/makefile
                 src/se_yield/makefile
                 src/si/makefile
                 src/ucmh/makefile
                 src/user_replaceables/makefile
                 src_models/common/aux/makefile
                 src_models/common/sm/makefile
                 src_models/empty/ac/makefile
                 src_models/empty/an/makefile
                 src_models/empty/cr/makefile
                 src_models/empty/dd/makefile
                 src_models/empty/ge/makefile
                 src_models/empty/ini/makefile
                 src_models/empty/kd/makefile
                 src_models/empty/rd/makefile
                 src_models/empty/se_yield/makefile
                 src_models/empty/si/makefile
                 src_models/empty/user_replaceables/makefile
                 src_models/generic_decayingDM/cr/makefile
                 src_models/generic_decayingDM/dec/makefile
                 src_models/generic_decayingDM/ge/makefile
                 src_models/generic_decayingDM/ini/makefile
                 src_models/generic_decayingDM/user_replaceables/makefile
                 src_models/generic_wimp/an/makefile
                 src_models/generic_wimp/cr/makefile
                 src_models/generic_wimp/dd/makefile
                 src_models/generic_wimp/ge/makefile
                 src_models/generic_wimp/ini/makefile
                 src_models/generic_wimp/kd/makefile
                 src_models/generic_wimp/rd/makefile
                 src_models/generic_wimp/se_yield/makefile
                 src_models/generic_wimp/user_replaceables/makefile
                 src_models/mssm/ac/makefile
                 src_models/mssm/ac_bsg/makefile
                 src_models/mssm/an/makefile
                 src_models/mssm/an_1l/makefile
                 src_models/mssm/an_ib/makefile
                 src_models/mssm/an_ib2/makefile
                 src_models/mssm/an_ib3/makefile
                 src_models/mssm/an_sf/makefile
                 src_models/mssm/an_stu/makefile
                 src_models/mssm/an_yield/makefile
                 src_models/mssm/an_yield_casc/makefile
                 src_models/mssm/cr/makefile
                 src_models/mssm/dd/makefile
                 src_models/mssm/examples/makefile
                 src_models/mssm/examples_aux/makefile
                 src_models/mssm/ge/makefile
                 src_models/mssm/ge_cmssm/makefile
                 src_models/mssm/ge_slha/makefile
                 src_models/mssm/ini/makefile
                 src_models/mssm/kd/makefile
                 src_models/mssm/rd/makefile
                 src_models/mssm/se_yield/makefile
                 src_models/mssm/se_yield_casc/makefile
                 src_models/mssm/user_replaceables/makefile
                 src_models/mssm/xfeynhiggs/makefile
                 src_models/mssm/xhiggsbounds/makefile
                 src_models/mssm/xsuperiso/makefile
                 src_models/silveira_zee/an/makefile
                 src_models/silveira_zee/cr/makefile
                 src_models/silveira_zee/dd/makefile
                 src_models/silveira_zee/ge/makefile
                 src_models/silveira_zee/ini/makefile
                 src_models/silveira_zee/kd/makefile
                 src_models/silveira_zee/rd/makefile
                 src_models/silveira_zee/se_yield/makefile
                 src_models/silveira_zee/user_replaceables/makefile
                 src_models/vdSIDM/an/makefile
                 src_models/vdSIDM/cr/makefile
                 src_models/vdSIDM/dd/makefile
                 src_models/vdSIDM/ge/makefile
                 src_models/vdSIDM/ini/makefile
                 src_models/vdSIDM/kd/makefile
                 src_models/vdSIDM/rd/makefile
                 src_models/vdSIDM/se_yield/makefile
                 src_models/vdSIDM/si/makefile
                 src_models/vdSIDM/user_replaceables/makefile
                 src/makefile
                 src_models/empty/makefile
                 src_models/generic_decayingDM/makefile
                 src_models/generic_wimp/makefile
                 src_models/mssm/makefile
                 src_models/silveira_zee/makefile
                 src_models/vdSIDM/makefile])

# Check if not_for_distribution exists, and if it does add its makefiles to
# configure
AC_CHECK_FILE(not_for_distribution/misc/makefile.in, AC_CONFIG_FILES([not_for_distribution/misc/makefile]))

AC_OUTPUT
