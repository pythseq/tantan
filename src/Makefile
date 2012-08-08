CXXFLAGS = -O3 -Wall -W -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -Wold-style-cast
# -Wconversion

CFLAGS = -Wall

COBJ = CA_code/lambda_calculator.o

all: tantan

tantan: *.cc *.hh Makefile $(COBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ *.cc $(COBJ)

$(COBJ): CA_code/*.c CA_code/*.h Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ CA_code/lambda_calculator.c

clean:
	rm -f tantan $(COBJ)
