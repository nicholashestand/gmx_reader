# simple makefile for order program

all: gmx_reader.o libgmx_reader.a

gmx_reader.o: gmx_reader.cpp gmx_reader.h
	g++ -c gmx_reader.cpp -lm -lxdrfile -std=c++11

libgmx_reader.a: gmx_reader.o
	ar rvs libgmx_reader.a gmx_reader.o

# note that gmx_reader must be linked first, then xdrfile
main_test.exe: main_test.cpp
	g++ main_test.cpp -o main_test.exe -lgmx_reader -lxdrfile 

install:
	mv libgmx_reader.a /usr/local/lib/
	cp gmx_reader.h /usr/local/include/

clean:
	rm -f gmx_reader.o main_test.exe libgmx_reader.a
