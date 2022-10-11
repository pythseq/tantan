CXXFLAGS = -msse4 -O3 -Wall -g

all: ../bin/tantan

../bin/tantan: *.cc *.hh version.hh Makefile
	mkdir -p ../bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ *.cc -lz

clean:
	rm -f ../bin/tantan

VERSION1 = git describe --dirty
VERSION2 = echo '$Format:%d$ ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e ../.git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@

FORCE:
