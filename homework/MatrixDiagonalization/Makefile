
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt eigenfunctions.png

eigenfunctions.png: out.txt
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "Energy";\
	set title "Eigenfunctions";\
	plot \
	"$<" index 1 using 1:2 with points color red title "Wavefunction numerical k=1" \
	,"$<" index 1 using 1:3 with points color blue title "Wavefunction numerical k=2" \
	,"$<" index 1 using 1:4 with points color green title "Wavefunction numerical k=3" \
	,"$<" index 1 using 1:5 with points color purple title "Wavefunction numerical k=4" \
	,"$<" index 1 using 1:6 with lines color red title "Wavefunction analytical k=1" \
	,"$<" index 1 using 1:7 with lines color blue title "Wavefunction analytical k=2" \
	,"$<" index 1 using 1:8 with lines color green title "Wavefunction analytical k=3" \
	,"$<" index 1 using 1:9 with lines color purple title "Wavefunction analytical k=4" \
	' | pyxplot

out.txt:main
	./$< > $@

main: main.c print.c jacobi.c

.PHONEY:clean
clean:
	$(RM) *txt main *png
