CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt data.txt test.png


out.txt: main
	./$< > $@

main: main.o HomeworkTools.c Eigenfinderdefault.o Eigenfinderguess.o EigenfinderGSL.o VectorTools.o

data.txt: plotter
	./$< > $@

plotter: plotter.o HomeworkTools.c Eigenfinderdefault.o VectorTools.o

test.png: data.txt Makefile
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set axis x;\
	set fontsize 0.4;\
	set grid;\
	set xlabel "x";\
	set ylabel "f(x)";\
	set title "Krakteristiske polynomoum med funde egenvaerdier";\
	plot \
	"$<" index 0 using 1:2 with points title "Egenvaerdier" \
	,"$<" index 1 using 1:2 with lines color blue title "Karakteristiske polynomium"\
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) out.txt data.txt plotter main *png *.o
