CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default: out.txt cos.png epidemic.png line.png Threebody.png


out.txt: main
	./$< > $@

main: main.c DriverAndStepper.c TestFunctions.c

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

Threebody.png: Threebody.txt
	echo '\
	set term png;\
	set output "$@";\
	set key top;\
	set fontsize 0.5;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Threebody system";\
	plot \
	 "$<" using 2:3 with lines title "Body 1"\
	,"$<" using 4:5 with lines title "Body 2"\
	,"$<" using 6:7 with lines title "Body 3"\
	' | pyxplot



.PHONEY:clean
clean:
	$(RM) *txt main *png
