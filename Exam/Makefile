CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt


out.txt: main
	./$< > $@

main: main.o HomeworkTools.c Eigenfinderdefault.o Eigenfinderguess.o EigenfinderGSL.o VectorTools.o





.PHONEY:clean
clean:
	$(RM) out.txt main *png *.o
