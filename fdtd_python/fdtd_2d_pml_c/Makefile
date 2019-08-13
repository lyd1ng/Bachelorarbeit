CC=gcc
CFLAGS=-Wno-deprecated-declarations
LDFLAGS=
LDLIBS=-lOpenCL -lm
OBJ=main.o fdtd_helper.o

main: $(OBJ)
	$(CC) -o main $(LDFLAGS) $(OBJ) $(LDLIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<


clean:
	rm -f main $(OBJ)

.PHONY:
	clean
	render
