CC=icc
CFLAGS=-O3 -Wall -Werror -lgsl -lgslcblas

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

OBJ=funct.o eccd.o
eccd: $(OBJ)
	icc -o $@ $^ $(CFLAGS)

OBJ=funct.o tresp.o
tresp: $(OBJ)
	icc -o $@ $^ $(CFLAGS)

clean:  
	rm -f eccd tresp *.o
