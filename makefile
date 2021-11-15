CC=g++

CFLAGS=-c -msse4.1

all: main

main: main.o 
	$(CC) $(LIBDIR) main.o -o main

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

clean:
	rm -f *o main