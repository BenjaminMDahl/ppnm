CFLAGS = -Wall -Werror -O1 # wall=all warnings werror=all warnings are errors so we leave file -01 optimzastion option
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY:default

default: out.main.txt myintegration.png

out.main.txt:main
	./$< > $@

myintegration.png: out.main.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "Upper integration limit";\
	set ylabel "integration value";\
	set title "int log(x)/sqrt(x)";\
	plot \
	 "$<" using 1:2 with line title "int log(x)/sqrt(x)"\
	' | tee tmp.gpi | pyxplot


.PHONEY:clean
clean:
	$(RM) out* main myintegration.png

