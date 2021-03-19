CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default:out.txt gsl.txt LinearSpline.png


out.txt: main
	./$< > $@

gsl.txt: gsl
	./$< > $@


LinearSpline.png: out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear Spline";\
	plot \
	"$<" index 0 using 1:2 with points pointtype 4 title "data" \
	,"$<" index 1 using 1:2 with lines color red title "linear-spline"\
	,"$<" index 1 using 1:3 with lines color red title "linear-spline-integral"\
	' | pyxplot
QudraticSpline.png: out.txt gsl.txt

CubicSpline.png: out.txt gsl.txt


.PHONEY:clean
clean:
	$(RM) *txt main gsl *png
