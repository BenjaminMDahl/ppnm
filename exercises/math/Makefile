CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

default: out.txt
	cat  out.txt

out.txt: math math2
	./math > out.txt
	./math2 >> out.txt

math: math.c
