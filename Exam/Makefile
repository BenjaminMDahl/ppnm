CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt test.png


out.txt: main
	./$< > $@

main: main.o HomeworkTools.c Eigenfinderdefault.o Eigenfinderguess.o EigenfinderGSL.o VectorTools.o

test.png: out.txt Makefile
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.4;\
	set xlabel "x";\
	set ylabel "f(x)";\
	set title "Krakteristiske polynomoum med funde egenvaerdier";\
	plot \
	"$<" index 1 using 1:2 with points title "Egenvaerdier" \
	,"$<" index 2 using 1:2 with lines color blue title "Karakteristiske polynomium"\
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) out.txt main *png *.o
