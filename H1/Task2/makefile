
CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm

HEADERS = initfcc.h alpotential.h utils.h equilibration.h
OBJECTS = initfcc.o alpotential.o utils.o equilibration.o MD_main.o
PROGRAM = run

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

