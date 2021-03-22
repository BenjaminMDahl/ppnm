CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default:out.main.txt out.quadratic.txt gsl.txt LinearSpline.png CubicSpline.png


out.main.txt: main
	./$< > $@

out.quadratic.txt: quadratic
	./$< > $@

gsl.txt: gsl
	./$< > $@


LinearSpline.png: out.main.txt Makefile gsl.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top left;\
	set fontsize 0.7;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear Spline";\
	plot \
	"$<" index 0 using 1:2 with points pointtype 3 title "data" \
	,"$<" index 1 using 1:2 with lines color red title "linear-spline"\
	,"$<" index 1 using 1:3 with lines color red title "linear-spline-integral"\
	,"gsl.txt" index 1 using 1:2 with lines linetype 7 color blue title "linear-GSL"\
	,"gsl.txt" index 2 using 1:2 with lines linetype 7 color blue title "linear-GSL-integral"\
	' | pyxplot

QudraticSpline.png: out.txt gsl.txt

CubicSpline.png:  out.main.txt Makefile gsl.txt
	echo '\
	set term png;\
	set output "$@";\
	set fontsize 0.65;\
	set key top left;\
	set keycolumns 2;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Cubic Spline";\
	plot \
	"$<" index 0 using 1:2 with points pointtype 3 title "data" \
	,"$<" index 3 using 1:2 with lines color red title "linear-spline"\
	,"$<" index 3 using 1:3 with lines color red title "linear-spline-derivative"\
	,"$<" index 3 using 1:4 with lines color red title "linear-spline-integral"\
	,"gsl.txt" index 1 using 1:3 with lines linetype 7 color blue title "linear-GSL"\
	,"gsl.txt" index 3 using 1:2 with lines linetype 7 color blue title "linear-GSL-derivative"\
	,"gsl.txt" index 2 using 1:3 with lines linetype 7 color blue title "linear-GSL-integral"\
	' | pyxplot


.PHONEY:clean
clean:
	$(RM) *txt main gsl *png quadratic
