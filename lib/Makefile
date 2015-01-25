COMPILE := gcc -std=c99 -Wall -Wextra -pedantic -O3 -pipe -MMD -I../include -c
LIBRARY := libmandelbrot-numerics.a
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))

all: $(LIBRARY) pkgconfig/mandelbrot-numerics.pc

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(LIBRARY) pkgconfig/mandelbrot-numerics.pc

$(LIBRARY): $(OBJECTS)
	@echo "AR      $@" ; ar -rs $@ $^ || ( echo "ERROR   ar -rs $@ $^" && false )

%.o: %.c
	@echo "C       $<" ; $(COMPILE) -o $@ $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

pkgconfig/mandelbrot-numerics.pc: pkgconfig/mandelbrot-numerics.pc.in
	@echo "PC      $<" ; ( echo "prefix=`cd .. ; pwd`" ; cat pkgconfig/mandelbrot-numerics.pc.in ) > pkgconfig/mandelbrot-numerics.pc || ( echo 'ERROR    ( echo "prefix=`cd .. ; pwd`" ; cat pkgconfig/mandelbrot-numerics.pc.in ) > pkgconfig/mandelbrot-numerics.pc' && false )

.SUFFIXES:
.PHONY: all clean

-include $(DEPENDS)
