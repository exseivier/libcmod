
FLAGS=-Wall -g -O3 -fPIC -shared
LIB=lib

build: libdnads.so

libdnads.so: $(LIB)/dnads.h $(LIB)/dnads.c
	gcc $(FLAGS) -o libdnads.so $(LIB)/dnads.c

clear:
	@rm libdnads.so
