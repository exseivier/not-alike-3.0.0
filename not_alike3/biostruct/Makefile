SRC:= $(wildcard lib/*.c)
OBJ:= $(SRC:.c=.o)

build: $(OBJ)
	gcc -Wall -g -o test $(OBJ)
	gcc -Wall -g -fPIC -shared -o libdnah.so -lm $(OBJ)

%.o: %.c
	gcc -Wall -g -c -o $@ $<

clean:
	rm lib/*.o *.so test
