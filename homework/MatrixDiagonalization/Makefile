
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt eigenfunctions.png

eigenfunctions.png: out.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "Energy";\
	set title "Eigenfunctions";\
	plot \
	"$<" index 1 using 1:2 with lines color red title "sin(x)" \
	,"$<" index 1 using 1:3 with lines color blue title "2*sin(2*x)" \
	' | pyxplot

out.txt:main
	./$< > $@

.PHONEY:clean
clean:
	$(RM) *txt main *png
