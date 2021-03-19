CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default:out.txt gsl.txt


out.txt: main
	./$< > $@

gsl.txt: gsl
	./$< > $@


LinearSpline.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "myerf"\
	,"wikidata.txt" using 1:2 with points pointtype 2 title "tabulated data"\
	' | | pyxplot

QudraticSpline.png: out.txt gsl.txt

CubicSpline.png: out.txt gsl.txt


.PHONEY:clean
clean:
	$(RM) *txt main gsl *png
