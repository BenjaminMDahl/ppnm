CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt


out.txt: main
	./$< > $@

main: main.c solver.c Eigenfinderdefault.c Eigenfinderguess.c




.PHONEY:clean
clean:
	$(RM) *txt main *png exact
