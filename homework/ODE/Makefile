CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt cos.png epidemic.png line.png


out.txt: main
	./$< > $@

line.png: line.txt
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


cos.png: cos.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "t";\
	set ylabel "y";\
	set title "ODE of cos(t)";\
	plot \
	"$<" using 1:2 with lines color red title "y from ODE" \
	,"$<" using 1:3 with lines color blue title "dy/dt from ODE" \
	' | pyxplot

epidemic.png: Epidemic.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "Days";\
	set ylabel "Number of people";\
	set title "Epidemic";\
	plot \
	"$<" using 1:2 with lines color red title "ODE for S" \
	,"$<" using 1:3 with lines color green title "ODE for I" \
	,"$<" using 1:4 with lines color blue title "ODE for R" \
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) *txt main *png
