CXXFLAGS = -O3 -Wall -W -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -Wold-style-cast
# -Wconversion

CFLAGS = -Wall

COBJ = CA_code/lambda_calculator.o

tantan: *.cc *.hh Makefile $(COBJ)
	$(CXX) $(CXXFLAGS) -o $@ *.cc $(COBJ)

$(COBJ): CA_code/*.c CA_code/*.h Makefile
	$(CC) $(CFLAGS) -c -o $@ CA_code/lambda_calculator.c
