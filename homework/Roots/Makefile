CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt kvant.png


out.txt: main
	./$< > $@

exact.txt: exact
	./$< > $@

main: main.o linsolve.o newton.o testfunctions.o ODE.o



kvant.png: kvant.txt exact.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "r";\
	set ylabel "f(r)";\
	set title "Hydrogen";\
	plot \
	"$<" using 1:2 with lines color blue title "line from ODE"\
	,"exact.txt" using 1:2 with lines linetype 7 color red title "exact line"\
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) *txt main *png exact
