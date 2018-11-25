
CFLAGS ?= -g -O2
MPICC ?= mpicc

all: mr_th_nb.c
	$(MPICC) $(CFLAGS) -o mr_th_nb mr_th_nb.c

clean:
	rm -f mr_th_nb