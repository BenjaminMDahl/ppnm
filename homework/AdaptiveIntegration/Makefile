CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt


out.txt: main
	./$< > $@

main: main.o TestFunctions.o integrater.o gsl_fun.o





.PHONEY:clean
clean:
	$(RM) *txt main *png
