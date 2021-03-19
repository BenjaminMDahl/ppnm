CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default:out.txt gsl.txt


out.txt: main
	./$< > $@

gsl.txt: gsl
	./$< > $@

.PHONEY:clean
clean:
	$(RM) *txt main gsl
