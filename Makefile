CC=gcc
CFLAGS=-Wall -Werror -g -O3
LDFLAGS=-lm -lpthread

.PHONY: all clean

all: pso

pso: pso.c pso.h problems.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ pso.c

clean:
	rm -rf pso
