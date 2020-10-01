# --- variables defined by configure ---

SRC = ./src
PREFIX = ds
LIBDIRSUFFIX = 64

EXE = 
DEF = -D
NOUNDERSCORE = 0
GM2FULL = 0
SLHAPARA = 0
VT100 = 1

FC = gfortran
FFLAGS = -O3 -ffixed-line-length-none $(DEF)U77EXT=0

CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -ffast-math -m64 -g -DBIGENDIAN=0

CXX = g++
CXXFLAGS = -m64 -stdlib=libstdc++

QUAD = 
HASQUAD = 1
ML = 0
MCC = mcc
MCFLAGS = -st -b64
MCLIBS = -lpthread

LDFLAGS = -L/usr/lib/gcc/x86_64-linux-gnu/9/liblto_plugin.so -L/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/9/../../.. -lgfortran -lm -lgcc_s -lgcc -lquadmath -lm -lgcc_s -lgcc -lgcc_s -lgcc -m64

AR = ar
RANLIB = ranlib
NM = nm
DLLTOOL = 

PARALLEL = -j 4

BLD = build

# --- end defs by configure ---


FHVERSION = 2.13.0

LIBDIR = $(PREFIX)/lib$(LIBDIRSUFFIX)
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

VERSION = $(BLD)/version.h

ARGS = $(PARALLEL) \
  GM2FULL="$(GM2FULL)" \
  NOUNDERSCORE="$(NOUNDERSCORE)" \
  EXE="$(EXE)" \
  FC="$(FC)" \
  FFLAGS="$(FFLAGS) $(DEF)HASQUAD=$(HASQUAD) $(DEF)GM2FULL=$(GM2FULL) $(DEF)SLHAPARA=$(SLHAPARA) $(DEF)VT100=$(VT100)" \
  CC="$(CC)" \
  CFLAGS="$(CFLAGS) -DNOUNDERSCORE=$(NOUNDERSCORE) -DVT100=$(VT100)" \
  CXX="$(CXX)" \
  CXXFLAGS="$(CXXFLAGS)" \
  QUAD="$(QUAD)" \
  ML="$(ML)" \
  MCC="$(MCC)" \
  MCFLAGS="$(MCFLAGS)" \
  MCLIBS="$(MCLIBS)" \
  LDFLAGS="$(LDFLAGS)" \
  AR="$(AR)" \
  RANLIB="$(RANLIB)" \
  NM="$(NM)" \
  DLLTOOL="$(DLLTOOL)"


default all lib frontend mma: force
	cd $(BLD) && $(MAKE) $(ARGS) $@

lib_shared: lib
	$(FC) $(FFLAGS) -shared -o $(BLD)/libFH.so build/*.o

install: default
	-mkdir "$(PREFIX)"
	-mkdir "$(LIBDIR)" "$(INCLUDEDIR)" "$(BINDIR)"
	cp -p $(BLD)/libFH.a "$(LIBDIR)"
	cp -p $(BLD)/CFeynHiggs.h $(BLD)/CSLHA.h $(BLD)/FHCouplings.h \
	  $(BLD)/FHRecord.h $(BLD)/SLHA.h $(BLD)/SLHADefs.h \
	  $(BLD)/PDG.h "$(INCLUDEDIR)"
	strip $(BLD)/FeynHiggs$(EXE) $(BLD)/table$(EXE)
	cp -p $(BLD)/FeynHiggs$(EXE) $(BLD)/table$(EXE) $(BLD)/fcc "$(BINDIR)"
	-strip $(BLD)/MFeynHiggs$(EXE)
	-cp -p $(BLD)/MFeynHiggs$(EXE) "$(BINDIR)"

force: $(BLD)/version.h

$(BLD)/version.h:
	-mkdir $(BLD)
	find "`cd $(SRC) ; pwd`" -follow -exec ln -sf {} $(BLD) \;
	( echo "#define FHVERSION \"$(FHVERSION)\"" ; \
	  echo "#define FHBUILDDATE \"`date '+%b %d, %Y'`\"" ) > $(BLD)/version.h


FH = FeynHiggs-$(FHVERSION)
FHTAR = $(FH).tar.gz

dist:
	rm -f $(FHTAR)
	ln -s . $(FH)
	tar cvfz $(FHTAR) --owner=root --group=root \
	  $(FH)/README $(FH)/configure $(FH)/makefile.in \
	  $(FH)/src $(FH)/extse $(FH)/gen \
	  $(FH)/man $(FH)/example $(FH)/web \
	  $(FH)/COPYING
	rm $(FH)

clean:
	rm -fr $(BLD) $(FHTAR)

.PHONY: default all install force clean dist

