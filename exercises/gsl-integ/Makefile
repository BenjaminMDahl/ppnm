CFLAGS = -Wall -Werror -O1 # wall=all warnings werror=all warnings are errors so we leave file -01 optimzastion option
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY:default

default: data.txt myintegration.png

data.txt:main
	./$< > $@

myintegration.png: data.txt Makefile Wolf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set fontsize 0.5;\
	set tics out;\
	set xlabel "Upper integration limit";\
	set ylabel "integration value";\
	set title "int log(x)/sqrt(x)";\
	plot \
	 "$<" using 1:2 with line title "GSL integral log(x)/sqrt(x)"\
	,"Wolf.txt" using 1:2 with points pointtype 3 title "Calculations from Wolframalpha.com"\
	' | tee tmp.gpi | pyxplot


.PHONEY:clean
clean:
	$(RM) out* data* main myintegration.png tmp.gpi

