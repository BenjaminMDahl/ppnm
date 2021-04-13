CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default:out.main.txt out.quadratic.txt gsl.txt LinearSpline.png CubicSpline.png QuadraticSpline.png


out.main.txt: main
	./$< > $@

main: main.c


out.quadratic.txt: quadratic
	./$< > $@

gsl.txt: gsl
	./$< > $@


LinearSpline.png: out.main.txt Makefile gsl.txt
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear Spline";\
	plot \
	"$<" index 0 using 1:2 with points pointtype 3 title "data(sin(x))" \
	,"$<" index 1 using 1:2 with lines color black title "linear-spline"\
	,"$<" index 1 using 1:3 with lines color green title "linear-spline-integral"\
	,"gsl.txt" index 1 using 1:2 with lines linetype 7 color blue title "linear-GSL"\
	,"gsl.txt" index 2 using 1:2 with lines linetype 7 color blue title "linear-GSL-integral"\
	' | pyxplot

QuadraticSpline.png: out.main.txt gsl.txt
	echo '\
	set term png;\
	set output "$@";\
	set fontsize 0.5;\
	set key outside;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Quadratic Spline";\
	plot \
	"$<" index 0 using 1:2 with points pointtype 3 title "data(sin(x))" \
	,"$<" index 2 using 1:2 with lines color black title "quadratic-spline"\
	,"$<" index 2 using 1:3 with lines color red title "quadratic-spline-diff"\
	,"$<" index 2 using 1:4 with lines color green title "quadratic-spline-int"\
	,"$<" index 0 using 1:3 with lines linetype 7 color blue title "analytic diff(cos(x))"\
	,"$<" index 0 using 1:4 with lines linetype 7 color blue title "analytic-int(-cos(x)+(k=1))"\
	' | pyxplot


CubicSpline.png:  out.main.txt Makefile gsl.txt
	echo '\
	set term png;\
	set output "$@";\
	set fontsize 0.5;\
	set key outside;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Cubic Spline";\
	plot \
	"$<" index 0 using 1:2 with points pointtype 3 title "data(sin(x))" \
	,"$<" index 3 using 1:2 with lines color black title "cubic-spline"\
	,"$<" index 3 using 1:3 with lines color red title "cubic-spline-diff"\
	,"$<" index 3 using 1:4 with lines color green title "cubic-spline-int"\
	,"gsl.txt" index 1 using 1:3 with lines linetype 7 color blue title "cubic-GSL"\
	,"gsl.txt" index 3 using 1:2 with lines linetype 7 color blue title "cubic-GSL-diff"\
	,"gsl.txt" index 2 using 1:3 with lines linetype 7 color blue title "linear-GSL-int"\
	' | pyxplot


.PHONEY:clean
clean:
	$(RM) *txt main gsl *png quadratic
