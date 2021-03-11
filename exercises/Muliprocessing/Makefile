CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
LDLIBS = -lm -lgomp



figure.png: points.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Monte Carlo";\
	plot \
	"$<" using 3:4 with line title "random numbers"\
	,"$<" using 1:2 with points pointtype 2 title "qudrant"\
	' | tee tmp.gpi | pyxplot


points.txt: numbers
	./$< > $@

numbers: numbers.o


.PHONEY:clean
clean:
	$(RM) *.o *.png *.txt numbers tmp.gpi

