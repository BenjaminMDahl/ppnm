CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

out.txt:main
	./$< > out.txt

main: main.o solver.o

.PHONEY:clean
clean:
	$(RM) *txt main
