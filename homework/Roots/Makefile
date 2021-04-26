CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt


out.txt: main
	./$< > $@

main: main.c linsolve.c newton.c testfunctions.c

line.png: out.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "t";\
	set ylabel "y";\
	set title "ODE of y(t)=t+1";\
	plot \
	"$<" using 1:2 with lines color blue title "line from ODE" \
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) *txt main *png
