CC = gcc # the C compiler
CFLAGS = -O -std=gnu11 # options for the C compiler
LDLIBS = -lm # libraries to link

default: out.txt
	cat out.txt

out.txt: hello
	./hello > out.txt

.PHONEY: clean
clean:                              # clean is a phoney target
	$(RM) hello.o hello out.txt # clean up the directory
