COMPILE := gcc -std=c99 -Wall -Wextra -pedantic -fPIC -O3 -pipe -MMD `PKG_CONFIG_PATH=../lib/pkgconfig pkg-config --cflags mandelbrot-numerics`
LIBS    := `PKG_CONFIG_PATH=../lib/pkgconfig pkg-config --libs mandelbrot-numerics`
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))
EXES    := $(patsubst %.o,%,$(OBJECTS))

all: $(EXES)

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(EXES)

%: %.o
	@echo "EXE     $@" ; gcc -o $@ $< $(LIBS) || ( echo "ERROR   gcc -o $@ $< $(LIBS)" && false )

%.o: %.c
	@echo "O       $@" ; $(COMPILE) -o $@ -c $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

.SUFFIXES:
.PHONY: all clean
.SECONDARY: $(OBJECTS)

-include $(DEPENDS)
