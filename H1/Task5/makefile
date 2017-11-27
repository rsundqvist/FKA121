
CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm

HEADERS = initfcc.h alpotential.h equilibration.h stat.h 
OBJECTS = initfcc.o alpotential.o equilibration.o stat.o MD_main.o
PROGRAM = run

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

