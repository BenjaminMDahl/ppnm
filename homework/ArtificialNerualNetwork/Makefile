CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt fit.png

out.txt: main
	./$< > $@

fit.png: out.txt Makefile
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "f(x)";\
	set title "Fit by neural network fitting";\
	plot \
	"$<" index 0 using 1:2 with points color red title "data points from function" \
	,"$<" index 1 using 1:2 with lines color blue title "line from ANN" \
	' | pyxplot

main: main.o qnewton.o network.o


clean:
	rm -f main *.o out* *png
