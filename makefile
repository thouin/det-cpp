.PHONY: all clean

CXXFLAGS:=-fopenmp -Wall -Wextra -O2
LIBS:=
#LIBS += -lgmp -lgmpxx

all: det

%.o: %.cpp det.hpp
	g++ $(CXXFLAGS) -c $< -o $@

det: main.o matrix.o rings.o exterior.o
	g++ $(CXXFLAGS) $^ $(LIBS) -o det

clean:
	rm -f det *.o
