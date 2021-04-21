CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

default: out.txt
	cat  out.txt

out.txt: math
	./math > out.txt


.PHONEY: clean
clean:
	$(RM) *.txt math
