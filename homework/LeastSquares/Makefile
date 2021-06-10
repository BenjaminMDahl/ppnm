CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt fit.png

fit.png: out.txt Makefile
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Least square";\
	plot \
	"$<" index 1 using 1:3:5 with errorbars title "data with log(y)" \
	,"$<" index 2 using 1:2 with lines color black title "least square line"\
	,"$<" index 2 using 1:3 with lines color red title "+ sigma"\
	,"$<" index 2 using 1:4 with lines color red title "- sigma"\
	,"$<" index 2 using 1:5 with lines color blue title "+ 2sigma"\
	,"$<" index 2 using 1:6 with lines color blue title "- 2sigma"\
	' | pyxplot

out.txt:main
	./$< > out.txt

main: main.o print.o solver.o LeastSquare.o

.PHONEY:clean
clean:
	$(RM) *txt main *png *.o
