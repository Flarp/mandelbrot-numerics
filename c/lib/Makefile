COMPILE := emcc -std=c99 -fPIC -O3 -pipe -MMD -I../include -I/usr/local/include -c
LINK    := emcc -shared 
LIBRARY := libmandelbrot-numerics
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))

all: $(LIBRARY).a $(LIBRARY).so pkgconfig/mandelbrot-numerics.pc

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(LIBRARY).a $(LIBRARY).so pkgconfig/mandelbrot-numerics.pc

install: $(LIBRARY).a $(LIBRARY).so ../include/mandelbrot-numerics.h pkgconfig/mandelbrot-numerics.pc
	install -d "$(prefix)/include" "$(prefix)/lib" "$(prefix)/lib/pkgconfig"
	install -m 644 -t "$(prefix)/include" ../include/mandelbrot-numerics.h
	install -m 644 -t "$(prefix)/lib" $(LIBRARY).a $(LIBRARY).so
	install -m 644 -t "$(prefix)/lib/pkgconfig" pkgconfig/mandelbrot-numerics.pc

$(LIBRARY).a: $(OBJECTS)
	@echo "A       $@" ; ar -rs $@ $^ || ( echo "ERROR   ar -rs $@ $^" && false )

$(LIBRARY).so: $(OBJECTS)
	@echo "SO      $@" ; $(LINK) -o $@ $^ -lmpc -lmpfr -lgmp -lm || ( echo "ERROR   $(LINK) -o $@ $^ -lmpc -lmpfr -lgmp -lm" && false )

%.o: %.c
	@echo "O       $@" ; $(COMPILE) -o $@ $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

pkgconfig/mandelbrot-numerics.pc: pkgconfig/mandelbrot-numerics.pc.in
	@echo "PC      $@" ; ( echo "prefix=$(prefix)" ; cat pkgconfig/mandelbrot-numerics.pc.in ) > pkgconfig/mandelbrot-numerics.pc || ( echo 'ERROR    ( echo "prefix=$(prefix)" ; cat pkgconfig/mandelbrot-numerics.pc.in ) > pkgconfig/mandelbrot-numerics.pc' && false )

.SUFFIXES:
.PHONY: all clean install

-include $(DEPENDS)
