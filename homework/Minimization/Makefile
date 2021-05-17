CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt BreitWigner.png

out.txt: main
	./$< > $@

BreitWigner.png: plot.txt data.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "Energy(GeV/(c*c))";\
	set ylabel "Cross Section";\
	set title "Fit to BreitWigner function";\
	plot \
	"plot.txt" using 1:2 with lines color blue title "line from fit" \
	,"data.txt" using 1:2 with points color red title "data from CERN" \
	' | pyxplot

main: main.c newton.c testfunctions.c

.PHONEY:clean
clean:
	$(RM) out.txt plot.txt main *png
