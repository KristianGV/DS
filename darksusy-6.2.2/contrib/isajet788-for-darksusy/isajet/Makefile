SOURCES := $(wildcard *.f)
OBJECTS = $(SOURCES:.f=.o)

# FFLAGS is set in the top-level Makefile
# -I directories are separate to get relative path correct
IFLAGS = -I../isainc

all: libisajet.a

libisajet.a: $(OBJECTS)
	ar r libisajet.a $(OBJECTS)

aldata.o: aldata.f
	$(FC) -c $(FFLAGS) $(IFLAGS) aldata.f

.f.o:
	$(FC) -c $(FFLAGS) $(IFLAGS) $<

clean:
	rm *.o libisajet.a
