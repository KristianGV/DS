SOURCES := $(wildcard *.tex)

all: isajet.ps

isajet.ps: isajet.dvi
	dvips isajet -o

isajet.dvi: $(SOURCES)
	latex isajet
	latex isajet
	latex isajet

clean:
	rm isajet.dvi isajet.ps isajet.aux isajet.log

