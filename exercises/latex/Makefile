CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exp.pdf: exp.txt exp-pyxplot.pdf
	pdflatex $<
	pdflatex $<

out.txt: main
	./$< > out.txt

exp-pyxplot.pdf: out.txt Makefile
	echo '\
	set term pdf;\
	set width 6;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "exp(x)";\
	set title "Exponential function";\
	plot \
	 "$<" using 1:2 with line title "exp from math.h"\
	,"$<" using 1:3 with points pointtype 2 title "ex from exercise"\
	' | pyxplot

.PHONEY: clean
clean:
	$(RM) exp.pdf exp.log exp.aux main out.txt exp.tex exp-pyxplot.pdf
