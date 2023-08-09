FLAGS=-Wall -g -fPIC -shared
LIB=lib
SRC=dnah.c
OUT=libdnah.so

build: $(OUT)

$(OUT): $(LIB)/$(SRC)
	gcc $(FLAGS) -o $(OUT) $(LIB)/$(SRC) -lm

delete:
	@rm $(OUT)
