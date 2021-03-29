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
	"$<" index 0 using 1:3:5 with errorbars title "data with log(y)" \
	,"$<" index 1 using 1:2 with lines color blue title "least square line"\
	' | pyxplot

out.txt:main
	./$< > out.txt

.PHONEY:clean
clean:
	$(RM) *txt main *png
