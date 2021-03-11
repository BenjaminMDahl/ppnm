CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
LDLIBS = -lm -lgomp



figure.png: points.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Monte Carlo";\
	plot \
	"$<" using 1:2 with line title "quadrant"\
	,"$<" using 3:4 with points pointtype 1 title "random numbers"\
	' | tee tmp.gpi | pyxplot


points.txt: numbers
	./$< > $@

numbers: numbers.o
