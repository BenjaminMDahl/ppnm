CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt cos.png epidemic.png


out.txt: main
	./$< > $@

cos.png: path.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "Energy";\
	set title "Eigenfunctions";\
	plot \
	"$<" using 1:2 with lines color red title "ODE of u''=-u" \
	' | pyxplot

epidemic.png: Epidemic.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "Energy";\
	set title "Eigenfunctions";\
	plot \
	"$<" using 1:2 with lines color red title "ODE for epidemic" \
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) *txt main *png
