CXXFLAGS = -O3 -Wall -W -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -Wold-style-cast
# -Wconversion

all: tantan

tantan: *.cc *.hh version.hh Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ *.cc

clean:
	rm -f tantan

VERSION = \"`hg id -n`\"

version.hh: FORCE
	if test -e ../.hg; \
	then echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@ ; \
	fi

FORCE:
