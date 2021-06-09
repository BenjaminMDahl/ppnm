CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: gamma.pyxplot.png


gamma.pyxplot.png: out.data.txt Makefile factor.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "n";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "gamma from math.h"\
	,"$<" using 1:3 with points pointtype 1 color black title "gsl-sf-gamma"\
	,"$<" using 1:4 with points pointtype 2 title "mygamma"\
	,"factor.data.txt" using 1:2 with points pointtype 3 color red title "tabulated data"\
	' | tee tmp.gpi | pyxplot

out.data.txt: main
	./$< > $@

main: main.o mygamma.o
main.o mygamma.o: mygamma.h

.PHONEY:clean
clean:
	$(RM) *.o main *.png out.data.txt tmp.gpi

