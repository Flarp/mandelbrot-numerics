COMPILE := gcc -std=c99 -Wall -Wextra -pedantic -fPIC -O3 -pipe -MMD -I../include -c
LINK    := gcc -shared
LIBRARY := libmandelbrot-numerics
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))

all: $(LIBRARY).a $(LIBRARY).so pkgconfig/mandelbrot-numerics.pc

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(LIBRARY).a $(LIBRARY).so pkgconfig/mandelbrot-numerics.pc

$(LIBRARY).a: $(OBJECTS)
	@echo "A       $@" ; ar -rs $@ $^ || ( echo "ERROR   ar -rs $@ $^" && false )

$(LIBRARY).so: $(OBJECTS)
	@echo "SO      $@" ; $(LINK) -o $@ $^ -lmpc -lmpfr -lgmp -lm || ( echo "ERROR   $(LINK) -o $@ $< -lmpc -lmpfr -lgmp -lm" && false )

%.o: %.c
	@echo "O       $@" ; $(COMPILE) -o $@ $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

pkgconfig/mandelbrot-numerics.pc: pkgconfig/mandelbrot-numerics.pc.in
	@echo "PC      $@" ; ( echo "prefix=`cd .. ; pwd`" ; cat pkgconfig/mandelbrot-numerics.pc.in ) > pkgconfig/mandelbrot-numerics.pc || ( echo 'ERROR    ( echo "prefix=`cd .. ; pwd`" ; cat pkgconfig/mandelbrot-numerics.pc.in ) > pkgconfig/mandelbrot-numerics.pc' && false )

.SUFFIXES:
.PHONY: all clean

-include $(DEPENDS)
