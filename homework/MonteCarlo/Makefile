CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt comparison.png

comparison.png: data.txt
	echo '\
	set term png;\
	set output "$@";\
	set key outside;\
	set fontsize 0.5;\
	set xlabel "N";\
	set ylabel "error";\
	set title "Quasi vs plain";\
	a=1; b=-0.5; f(x)=a*x**b ;\
	c=1; d=-1.0; g(x)=c*x**d ;\
	fit f(x) "$<" index 0 using 1:2 via a,b;\
	fit g(x) "$<" index 0 using 1:3 via c,d;\
	plot \
	"$<" index 0 using 1:2 with points pointtype 2 title "Error from plainMC" \
	,"$<" index 0 using 1:3 with points pointtype 3 title "Error from QuasiMC"\
	,f(x) with lines title "f(x)" \
	,g(x) with lines title "g(x)" \
	' | pyxplot


data.txt: cirkel
	./cirkel > data.txt

out.txt: main
	./main > out.txt

main: main.c montecarlo.c

cirkel: cirkel.c montecarlo.c

.PHONEY:clean
clean:
	$(RM) *txt main *png cirkel
