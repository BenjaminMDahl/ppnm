
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp


figure.png: points.txt Makefile cir.data
	echo '\
	set terminal png;\
	set output "$@";\
	set key outside;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Monte Carlo";\
	plot \
	"$<" using 1:2 with points pointtype 1 pointsize 0.1 title "Random points"\
	,"cir.data" using 1:2 with line title "Quadrand"\
	' | tee tmp.gpi | pyxplot
	cat out.txt

points.txt: numbers
	./$< > $@

numbers: numbers.o


.PHONEY:clean
clean:
	$(RM) *.o *.png *.txt numbers tmp.gpi

