CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY:default
default: erf.pyxplot.png


erf.pyxplot.png: Makefile wikidata.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "DATA";\
	plot \
	 "wikidata.txt" using 1:2 with points pointtype 2 title "tabulated data"\
	' | tee tmp.gpi | pyxplot

tmp.gpi:erf.pyxplot.png
	touch $@

.PHONEY:clean
clean:
	$(RM) *.o main *.png out.data.txt tmp.gpi
