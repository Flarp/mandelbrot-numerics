PKGCONFIG := PKG_CONFIG_PATH="../lib/pkgconfig" pkg-config
COMPILE := emcc -I../include -I/usr/local/include -std=c99 -Wall -Wextra -pedantic -fPIC -O3 -pipe -MMD `$(PKGCONFIG) --cflags mandelbrot-numerics`
LIBS    := `$(PKGCONFIG) --libs mandelbrot-numerics` -lmpc -lmpfr -lgmp -lm
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))
EXES    := $(patsubst %.o,%,$(OBJECTS))

all: $(EXES)

out: 
	@echo $(PKGCONFIG)

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(EXES)

install: $(EXES)
	install -d "$(prefix)/bin"
	install -m 755 -t "$(prefix)/bin" $(EXES)

%: %.o
	@echo "EXE     $@" ; emcc -o $@ $< $(LIBS) || ( echo "ERROR   emcc -o $@ $< $(LIBS)" && false )

%.o: %.c
	@echo "O       $@" ; $(COMPILE) -o $@ -c $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

.SUFFIXES:
.PHONY: all clean install
.SECONDARY: $(OBJECTS)

-include $(DEPENDS)
