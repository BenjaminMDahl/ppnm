CFLAGS = -o -w -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: output.txt ; cat $<
output.txt: main ; ./$< > $@
clean: ; $(RM) main output.txt
